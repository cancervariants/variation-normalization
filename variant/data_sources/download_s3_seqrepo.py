"""Elastic Beanstalk to download seqrepo from s3."""
import boto3
import os
import zipfile
import shutil

s3 = boto3.resource('s3')
zip_path = "temp_seqrepo.zip"
bucket = os.environ['AWS_BUCKET_NAME']
obj = os.environ['AWS_SEQREPO_OBJECT']
s3.meta.client.download_file(bucket, obj, zip_path)
data_dir = "../data"
with zipfile.ZipFile(zip_path, 'r') as zip_ref:
    zip_ref.extractall(data_dir)
os.remove(zip_path)
shutil.rmtree(f"{data_dir}/__MACOSX")
