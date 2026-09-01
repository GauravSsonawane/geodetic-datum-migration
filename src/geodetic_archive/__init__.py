"""Migrate a legacy geodetic archive onto WGS84.

The seven-parameter transformation between the two datums is not published, so it
is recovered from the archive's own connection stations before anything can be
migrated. See :mod:`geodetic_archive.pipeline` for the end-to-end run.
"""

__version__ = "0.1.0"
