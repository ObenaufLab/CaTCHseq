#!/usr/bin/env python3


""" 
Provide a command line tool to validate and transform tabular samplesheets.
See https://github.com/ObenaufLab/CaTCHseq/blob/dev/bin/check_samplesheet.py for template
"""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (".fastq.gz",)

    def __init__(
        self,
        sample_col="SampleName",
        cond_col="Condition",
        rep_col="Replicate",
        lib_col="LibraryType",
        R1_col="R1",
        R2_col="R2",
        cellnr_col="CellNumber",
        chem_col="Chemistry",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "SampleName").
            cond_col (str): The name of the column that contains the experiments condition (e.g. treatment or timepoint, default "Condition0").
            rep_col (str): The name of the column that contains the replicate number (default "Replicate").
            lib_col (str): The name of the column that contains a library identifier (e.g. 10X,ScaleBio default "LibraryType").
            R1_col (str): The name of the column that contains the first FASTQ file path (default "R1").
            R2_col (str): The name of the column that contains the second FASTQ file path (default "R2").
            cellnr_col (str): The name of the column that contains the estimated cell number (e.g. 1000, NA, 1000!, default "CellNumber").
            chem_col (str): The name of the column that contains the sequencing chemistry (e.g. 10X,ScaleBio default "Chemistry").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._cond_col = cond_col
        self._rep_col = rep_col
        self._lib_col = lib_col
        self._R1_col = R1_col
        self._R2_col = R2_col
        self._cellnr_col = cellnr_col
        self._chem_col = chem_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_condition(row)
        self._validate_replicate(row)
        self._validate_library(row)
        self._validate_R1(row)
        self._validate_R2(row)
        self._validate_cellnr(row)
        self._validate_chemistry(row)
        self._seen.add(
            (
                row[self._sample_col],
                row[self._cond_col],
                row[self._rep_col],
                row[self._lib_col],
                row[self._R1_col],
            )
        )
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_condition(self, row):
        """Assert that the condition name exists and convert spaces to underscores."""
        if len(row[self._cond_col]) <= 0:
            raise AssertionError("Condition input is required.")
        # Sanitize condition slightly.
        row[self._cond_col] = row[self._cond_col].replace(" ", "_")

    def _validate_replicate(self, row):
        """Assert that the replicate exists and is a positive integer."""
        if len(row[self._rep_col]) <= 0:
            raise AssertionError("Replicate input is required.")
        try:
            replicate = int(row[self._rep_col])
            if replicate <= 0:
                raise AssertionError("Replicate number must be a positive integer.")
        except ValueError:
            raise AssertionError("Replicate number must be a positive integer.")

    def _validate_library(self, row):
        """Assert that the library type exists and is a valid value."""
        if len(row[self._lib_col]) <= 0:
            raise AssertionError("Library type input is required.")
        # Add your validation logic for library type here.

    def _validate_R1(self, row):
        """Assert that the R1 FASTQ file path exists and has the right format."""
        if len(row[self._R1_col]) <= 0:
            raise AssertionError("R1 FASTQ file path is required.")
        # Add your validation logic for R1 FASTQ file path here.

    def _validate_R2(self, row):
        """Assert that the R2 FASTQ file path exists and has the right format if it exists."""
        if len(row[self._R2_col]) > 0:
            # Add your validation logic for R2 FASTQ file path here.
            pass

    def _validate_cellnr(self, row):
        """Assert that the estimated cell number exists and is a positive integer or NA."""
        if len(row[self._cellnr_col]) <= 0:
            raise AssertionError("Estimated cell number input is required.")
        # Add your validation logic for estimated cell number here.

    def _validate_chemistry(self, row):
        """Assert that the sequencing chemistry exists and is a valid value."""
        if len(row[self._chem_col]) <= 0:
            raise AssertionError("Sequencing chemistry input is required.")
        # Add your validation logic for sequencing chemistry here.

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._R1_col] and row[self._R2_col]:
            first_col_suffix = Path(row[self._R1_col]).suffixes[-2:]
            second_col_suffix = Path(row[self._R2_col]).suffixes[-2:]
            if first_col_suffix != second_col_suffix:
                raise AssertionError("FASTQ pairs must have the same file extensions.")
        else:
            raise AssertionError("Read pairs are expected.")

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FORMATS):
            raise AssertionError(
                f"The FASTQ file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FORMATS)}"
            )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename all samples to have a suffix of _T{n}, where n is the
        number of times the same sample exist, but with different FASTQ files, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and FASTQ must be unique.")
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample_col]
            seen[sample] += 1
            row[self._sample_col] = f"{sample}_T{seen[sample]}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    # Function does not work properly
    # if not sniffer.has_header(peek):
    #     logger.critical("The given sample sheet does not appear to contain a header.")
    #     sys.exit(1)
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in):
    """
    Check that the tabular samplesheet has the structure expected by our pipeline.

    Validate the general shape of the table, expected columns, and each row

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The expected format is CSV

    Example:
        This function checks that the samplesheet follows the following structure

            SampleName,Condition,Replicate,LibraryType,R1,R2,CellNumber,Chemistry
            ScaleBio,Scale,1,GEX,S1_R1_001.fastq.gz,S1_R2_001.fastq.gz,NA,ScaleBio
            ScaleBio,Scale,1,scCaTCH,S6_L001_R1_001.fastq.gz,S6_L001_R2_001.fastq.gz,NA,ScaleBio
            Day0,Day0,1,GEX,S2_L001_R1_001.fastq.gz,S2_L001_R2_001.fastq.gz,7500!,10X
            Day0,Day0,1,scCaTCH,S3_L001_R1_001.fastq.gz,S3_L001_R2_001.fastq.gz,NA,10X

    """
    required_columns = {
        "SampleName",
        "Condition",
        "Replicate",
        "LibraryType",
        "R1",
        "R2",
        "CellNumber",
        "Chemistry",
    }
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(
                f"The sample sheet **must** contain these column headers: {req_cols}."
            )
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with open("sanity_check", mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate a CSV samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    check_samplesheet(args.file_in)


if __name__ == "__main__":
    sys.exit(main())  #!/usr/bin/env python
