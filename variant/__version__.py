"""Module for version information."""
__title__ = 'variant normalization'
__description__ = 'The VICC Variant Normalizer'
__url__ = 'http://normalize.cancervariants.org/variant'
__major__ = '0'
__minor__ = '0'
__patch__ = '1'
__meta_label__ = ''
__short_version__ = "{}.{}".format(__major__, __minor__)
__version__ = "{}.{}".format(__short_version__, __patch__)
if __meta_label__:
    __version__ += "-{}".format(__meta_label__)
__authors__ = ['Alex H. Wagner', 'Adam Coffman']
__author_email__ = 'help@variant-normalization.org'  # TODO: Fix this
__license__ = 'MIT'
__copyright__ = 'Copyright 2020 The Wagner Lab'
