"""A module for accessing SeqRepo."""
from typing import Optional
from biocommons.seqrepo import SeqRepo
from variant import SEQREPO_DATA_PATH, PROJECT_ROOT
import os
import boto3
import zipfile
import shutil


class SeqRepoAccess:
    """The SeqRepoAccess class."""

    def __init__(self, seqrepo_data_path=SEQREPO_DATA_PATH):
        """Initialize the SeqRepoAccess class.

        :param str seqrepo_data_path: The path to the seqrepo directory.
        """
        if not os.path.exists:
            self._download_from_s3()
        self.seq_repo_client = SeqRepo(seqrepo_data_path)

    def _download_from_s3(self):
        """Download SeqRepo data for Elastic Beanstalk."""
        # TODO: Consider doing this in .ebextensions?
        s3 = boto3.resource('s3')
        zip_path = f"{PROJECT_ROOT}/test.zip"
        bucket = os.environ['AWS_BUCKET_NAME']
        obj = os.environ['AWS_SEQREPO_OBJECT']
        s3.meta.client.download_file(bucket, obj, zip_path)
        data_dir = f"{PROJECT_ROOT}/variant/data/seqrepo"
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(data_dir)
        os.remove(zip_path)
        shutil.rmtree(f"{data_dir}/__MACOSX")

    def protein_at_position(self, transcript: str, pos: int) -> Optional[str]:
        """Get the protein at a position."""
        # why does this not exist sometimes?
        try:
            t = self.seq_repo_client.fetch(transcript)
            if len(t) < pos - 1:
                return None
            else:
                try:
                    return t[pos - 1]
                except IndexError:
                    return None
        except KeyError:
            return None

    def aliases(self, input_str):
        """Get aliases for gene symbols."""
        try:
            return self.seq_repo_client.translate_alias(input_str.strip())
        except KeyError:
            return []
