#!/usr/bin/env python
"""prepare_hetsnps_for_frankenplot.py

Quality-filter a VCF, keep only heterozygous simple SNVs, and add a new sample column
carrying DP and AD (ref,alt observations) read from a BAM pileup at each site.

Usage example:
    prepare_hetsnps_for_frankenplot.py -i germline.vcf.gz -b tumor.bam -o out.vcf.gz \\
        --samplename TUMOR --filter-hom --no-filtered --site-quality 5

VCF and BAM parsing use pysam; the bgzipped output and its tabix index are written by this
script.
"""

import argparse
import logging
import sys
from collections import Counter
from collections.abc import Callable, Generator, Sequence
from contextlib import contextmanager

import pysam
from pysam import (
    AlignmentFile,
    BGZFile,
    VariantFile,
    VariantHeader,
    VariantRecord,
    VariantRecordSample,
)

__version__ = "1.0.0"

logger = logging.getLogger("prepare_hetsnps_for_frankenplot")

MISSING = "."
NOCALL_GT = "./."
DEFAULT_MAX_DEPTH = 8000
SITE_QUALITY_FILTER_DESCRIPTION = "Filter low quality sites"

# FORMAT fields written for the added sample; used as fallback header definitions when the
# input header does not already define them.
FORMAT_DEFINITIONS: dict[str, str] = {
    "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
    "AD": (
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description='
        '"Allelic depths for the ref and alt alleles in the order listed">'
    ),
}


def is_homozygous(sample: VariantRecordSample) -> bool:
    """True when every called allele is the same, e.g. 0/0, 1/1 or the phased 0|0."""
    indices = sample.allele_indices
    return bool(indices) and None not in indices and len(set(indices)) == 1


def site_quality_filter_id(site_quality: int) -> str:
    """FILTER id for a quality threshold, e.g. ``sq5`` (PyVCF's naming)."""
    return f"sq{site_quality}"


def is_low_quality(record: VariantRecord, site_quality: int) -> bool:
    """True when QUAL is below the threshold.

    A missing QUAL ('.') passes: PyVCF on Python 2 compared None < threshold and returned
    None, which its filter loop read as "not filtered".
    """
    return record.qual is not None and float(record.qual) < site_quality


def is_snv(record: VariantRecord) -> bool:
    """True for a single biallelic 1bp REF/ALT swap; indels and multi-allelics are not."""
    alts = record.alts or ()
    return len(record.ref or "") == 1 and len(alts) == 1 and len(alts[0]) == 1


def pileup_bases(bam: AlignmentFile, contig: str, position: int, max_depth: int) -> list[str]:
    """Collect the read bases aligned to a 1-based position, skipping deletions and skips.

    pysam's pileup defaults apply, as they did in vcf_add_sample.py: reads flagged unmapped,
    secondary, QC-fail or duplicate are skipped, bases below quality 13 are ignored, and
    overlapping mates are counted once. ``max_depth`` caps how many reads are piled up at a
    site; deep panels need it raised above the pysam default.
    """
    bases: list[str] = []
    for column in bam.pileup(contig, position - 1, position, truncate=True, max_depth=max_depth):
        for read in column.pileups:
            if read.is_del or read.is_refskip:
                continue
            sequence = read.alignment.query_sequence
            if sequence is None or read.query_position is None:
                continue
            bases.append(sequence[read.query_position])
    return bases


VcfTextWriter = Callable[[str], object]


@contextmanager
def open_output(output: str) -> Generator[VcfTextWriter, None, None]:
    """Yield a text writer for ``-`` (stdout), ``*.vcf`` (plain) or ``*.gz`` (BGZF)."""
    if output == "-":
        yield sys.stdout.write
        return
    if output.endswith(".gz"):
        bgzf = BGZFile(output, "wb")
        try:
            yield lambda text: bgzf.write(text.encode("utf-8"))
        finally:
            bgzf.close()
        return
    with open(output, "w", encoding="utf-8") as plain:
        yield plain.write


def build_output_header(template: VariantHeader, samplename: str) -> VariantHeader:
    """Copy the input header, declare the FORMAT keys written here, and add the sample."""
    header = template.copy()
    for key, definition in FORMAT_DEFINITIONS.items():
        if key not in header.formats:
            header.add_line(definition)
    header.add_sample(samplename)
    return header


def render_record(record: VariantRecord, new_sample: dict[str, str]) -> str:
    """Render a VCF line with an extra sample column appended.

    pysam cannot append a sample to a record parsed against an existing header, so the line
    is rendered by htslib and the new column is appended textually. FORMAT keys missing from
    the record (typically DP or AD) are added, and the existing samples padded with '.'.
    """
    fields = str(record).rstrip("\n").split("\t")
    fixed, rest = fields[:8], fields[8:]
    keys = rest[0].split(":") if rest else ["GT"]
    samples = rest[1:]

    for key in new_sample:
        if key not in keys:
            keys.append(key)
            samples = [f"{sample}:{MISSING}" for sample in samples]

    added = ":".join(new_sample.get(key, MISSING) for key in keys)
    return "\t".join([*fixed, ":".join(keys), *samples, added]) + "\n"


def process(
    vcf_path: str,
    bam_path: str,
    output: str,
    samplename: str,
    site_quality: int | None,
    filter_hom: bool,
    drop_filtered: bool,
    max_depth: int,
) -> int:
    """Filter, subset to het SNVs, add BAM-derived depths, and write the output VCF."""
    written = 0

    with VariantFile(vcf_path) as in_vcf, AlignmentFile(bam_path, "rb") as bam:
        if samplename in in_vcf.header.samples:
            raise ValueError(f"Sample '{samplename}' is already present in {vcf_path}")

        # The FILTER id has to be declared on the header the records are attached to:
        # htslib rejects record.filter.add() for an id its own header does not know.
        filter_id = site_quality_filter_id(site_quality) if site_quality is not None else None
        if filter_id is not None:
            in_vcf.header.filters.add(filter_id, None, None, SITE_QUALITY_FILTER_DESCRIPTION)
        header = build_output_header(in_vcf.header, samplename)

        with open_output(output) as writer:
            writer(str(header))
            for record in in_vcf:
                if site_quality is not None and is_low_quality(record, site_quality):
                    if drop_filtered:
                        continue
                    record.filter.add(filter_id)
                # PyVCF only stamped PASS when filtered records were kept in the output.
                # Pre-existing FILTER values from the input are left untouched either way.
                if not drop_filtered and not list(record.filter):
                    record.filter.add("PASS")

                if filter_hom and any(is_homozygous(call) for call in record.samples.values()):
                    continue

                if not is_snv(record):
                    continue

                alt = record.alts[0]
                counts = Counter(pileup_bases(bam, record.contig, record.pos, max_depth))
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(
                        "pileup at %s:%d %s/%s = %s",
                        record.contig,
                        record.pos,
                        record.ref,
                        alt,
                        "".join(base * count for base, count in sorted(counts.items())),
                    )
                writer(
                    render_record(
                        record,
                        {
                            "GT": NOCALL_GT,
                            "DP": str(sum(counts.values())),
                            "AD": f"{counts[record.ref]},{counts[alt]}",
                        },
                    )
                )
                written += 1

    logger.info("Wrote %d variant(s) to %s", written, output)
    if output.endswith(".gz"):
        pysam.tabix_index(output, preset="vcf", force=True)
        logger.info("Indexed %s.tbi", output)
    return written


def setup_logging(loglevel: str) -> None:
    """Configure stderr logging at the requested level.

    An unknown level name raises ValueError out of the logging module itself.
    """
    logging.basicConfig(
        level=loglevel.upper(),
        stream=sys.stderr,
        format="%(levelname)s %(asctime)s %(funcName)s - %(message)s",
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the command line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Filter a VCF, keep heterozygous simple SNVs, and add a sample column with "
            "DP and AD taken from a BAM pileup. Writes .vcf.gz plus its tabix index."
        ),
        epilog="Compression is single threaded (htslib BGZF).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Input VCF/BCF ('-' for stdin)")
    parser.add_argument("-b", "--bam", required=True, help="Indexed BAM for the sample to add")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output VCF; '.gz' writes BGZF and a tabix index, '-' writes to stdout",
    )
    parser.add_argument("--samplename", required=True, help="Name of the sample to add")
    parser.add_argument(
        "--site-quality",
        type=int,
        default=None,
        metavar="QUAL",
        help="Filter sites with QUAL below this value; omit to disable quality filtering",
    )
    parser.add_argument(
        "--filter-hom",
        action="store_true",
        help="Drop variants that are homozygous in any input sample",
    )
    parser.add_argument(
        "--no-filtered",
        action="store_true",
        help="Output only sites passing the filter, instead of tagging them in FILTER",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=DEFAULT_MAX_DEPTH,
        metavar="N",
        help="Maximum number of reads piled up at a site; raise it for deep panels",
    )
    parser.add_argument("--loglevel", default="INFO", help="Level of logging")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Parse arguments and run the conversion."""
    args = build_parser().parse_args(argv)
    setup_logging(args.loglevel)
    logger.info("Started log with loglevel %s", args.loglevel)

    process(
        vcf_path=args.input,
        bam_path=args.bam,
        output=args.output,
        samplename=args.samplename,
        site_quality=args.site_quality,
        filter_hom=args.filter_hom,
        drop_filtered=args.no_filtered,
        max_depth=args.max_depth,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
