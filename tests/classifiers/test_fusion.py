import unittest

from varlexapp.classifiers import FusionClassifier
from .classifier_base import ClassifierBase

class TestGenePairTokenizer(ClassifierBase, unittest.TestCase):

    def classifier_instance(self):
        return FusionClassifier()

    def fixture_name(self):
        return 'fusion'
