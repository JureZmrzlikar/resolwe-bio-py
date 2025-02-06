""".. Ignore pydocstyle D400.

========
QCTables
========

.. autoclass:: QCTables
    :members:
    :inherited-members:

    .. automethod:: __init__

"""

from functools import lru_cache
from typing import Callable, Dict, List, Optional

import pandas as pd

from resdk.resources import Data

from .base import BaseTables

CHUNK_SIZE = 1000


MQC_GENERAL_COLUMNS = [
    {
        "name": "FastQC (raw)_mqc-generalstats-fastqc_raw-total_sequences",
        "slug": "total_read_count_raw",
        "type": "Int64",
        "agg_func": "sum",
    },
    {
        "name": "FastQC (trimmed)_mqc-generalstats-fastqc_trimmed-total_sequences",
        "slug": "total_read_count_trimmed",
        "type": "Int64",
        "agg_func": "sum",
    },
    {
        "name": "FastQC (raw)_mqc-generalstats-fastqc_raw-percent_gc",
        "slug": "gc_content_raw",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "FastQC (trimmed)_mqc-generalstats-fastqc_trimmed-percent_gc",
        "slug": "gc_content_trimmed",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "FastQC (raw)_mqc-generalstats-fastqc_raw-percent_duplicates",
        "slug": "seq_duplication_raw",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "FastQC (trimmed)_mqc-generalstats-fastqc_trimmed-percent_duplicates",
        "slug": "seq_duplication_trimmed",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "FastQC (raw)_mqc-generalstats-fastqc_raw-avg_sequence_length",
        "slug": "avg_seq_length_raw",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "FastQC (trimmed)_mqc-generalstats-fastqc_trimmed-avg_sequence_length",
        "slug": "avg_seq_length_trimmed",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "STAR_mqc-generalstats-star-uniquely_mapped_percent",
        "slug": "mapped_reads_percent",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "STAR_mqc-generalstats-star-uniquely_mapped",
        "slug": "mapped_reads",
        "type": "Int64",
        "agg_func": "sum",
    },
    {
        "name": "STAR (Globin)_mqc-generalstats-star_globin-uniquely_mapped_percent",
        "slug": "mapped_reads_percent_globin",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "STAR (Globin)_mqc-generalstats-star_globin-uniquely_mapped",
        "slug": "mapped_reads_globin",
        "type": "Int64",
        "agg_func": "sum",
    },
    {
        "name": "STAR (rRNA)_mqc-generalstats-star_rrna-uniquely_mapped_percent",
        "slug": "mapped_reads_percent_rRNA",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "STAR (rRNA)_mqc-generalstats-star_rrna-uniquely_mapped",
        "slug": "mapped_reads_rRNA",
        "type": "Int64",
        "agg_func": "sum",
    },
    {
        "name": "featureCounts_mqc-generalstats-featurecounts-percent_assigned",
        "slug": "fc_assigned_reads_percent",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "featureCounts_mqc-generalstats-featurecounts-Assigned",
        "slug": "fc_assigned_reads",
        "type": "Int64",
        "agg_func": "sum",
    },
    {
        "name": "STAR quantification_mqc-generalstats-star_quantification-of_assigned_reads",
        "slug": "star_assigned_reads_percent",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "STAR quantification_mqc-generalstats-star_quantification-Assigned_reads",
        "slug": "star_assigned_reads",
        "type": "Int64",
        "agg_func": "sum",
    },
    {
        "name": "Salmon_mqc-generalstats-salmon-percent_mapped",
        "slug": "salmon_assigned_reads_percent",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "Salmon_mqc-generalstats-salmon-num_mapped",
        "slug": "salmon_assigned_reads",
        "type": "Int64",
        "agg_func": "sum",
    },
    {
        "name": "QoRTs_mqc-generalstats-qorts-Genes_PercentWithNonzeroCounts",
        "slug": "nonzero_count_features_percent",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "QoRTs_mqc-generalstats-qorts-NumberOfChromosomesCovered",
        "slug": "contigs_covered",
        "type": "Int64",
        "agg_func": "mean",
    },
    {
        "slug": "strandedness_code",
        "type": "string",
    },
    {
        "slug": "genome_build",
        "type": "string",
    },
]

QORTS_COLUMNS = [
    {
        "name": "StrandTest_frFirstStrand",
        "slug": "first_strand",
        "type": "float64",
        "agg_func": "mean",
    },
    {
        "name": "StrandTest_frSecondStrand",
        "slug": "second_strand",
        "type": "float64",
        "agg_func": "mean",
    },
]

def general_multiqc_parser(file_object, name, column_names):
    """General parser for MultiQC files."""
    df = pd.read_csv(file_object, sep="\t", index_col=0)

    # Keep only specified columns:
    df = df[
        [
            column.get("name", "")
            for column in column_names
            if column.get("name", "") in df.columns
        ]
    ]
    # Rename
    df = df.rename(
        columns={
            column.get("name", ""): column["slug"]
            for column in column_names
            if column.get("name", "") in df.columns
        }
    )

    if df.empty:
        return pd.Series(name=name)

    # Perform aggregation
    series = df.agg(
        {
            column["slug"]: column["agg_func"]
            for column in column_names
            if column["slug"] in df.columns
        }
    )
    series.name = name
    return series


class QCTables(BaseTables):
    """A helper class to fetch collection's QC data.

    A simple example:

    .. code-block:: python

        # Get Collection object
        collection = res.collection.get("collection-slug")

        # Fetch collection expressions and metadata
        tables = QCTables(collection)
        tables.qc

    """

    process_type = "data:multiqc:"

    # Data types:
    GENERAL = "general"
    QORTS = "qorts"

    data_type_to_field_name = {
        GENERAL: "report_data",
        QORTS: "report_data",
    }

    def _parse_file(self, file_obj, sample_id, data_type):
        """Parse file object and return a one DataFrame line."""
        if data_type == self.GENERAL:
            return general_multiqc_parser(file_obj, sample_id, MQC_GENERAL_COLUMNS)
        if data_type == self.QORTS:
            return general_multiqc_parser(file_obj, sample_id, QORTS_COLUMNS)


    def _get_data_uri(self, data: Data, data_type: str) -> str:
        if data_type == self.GENERAL:
            return f"{data.id}/multiqc_data/multiqc_general_stats.txt"
        if data_type == self.QORTS:
            return f"{data.id}/multiqc_data/multiqc_qorts.txt"

    @property
    @lru_cache()
    def general(self) -> pd.DataFrame:
        return self._load_fetch(self.GENERAL)

    @property
    @lru_cache()
    def qorts(self) -> pd.DataFrame:
        return self._load_fetch(self.QORTS)

    @property
    @lru_cache()
    def qc(self):
        # If needed, we can just concatenate all the outputs:
        return pd.concat(
            objs=[
                self.general,
                self.qorts,
                # ...
            ],
            axis=1,
        )
