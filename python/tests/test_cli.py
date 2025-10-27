"""Tests for CLI helper utilities."""

from bioscript.cli import _merge_classifier_result


def test_merge_classifier_result_string():
    results = {}
    _merge_classifier_result(results, "HERC2", "blue")

    assert results == {"HERC2_result": "blue"}


def test_merge_classifier_result_dict_with_result_key():
    results = {}
    _merge_classifier_result(results, "HERC2", {"result": "blue", "genotype": "GG"})

    assert results == {
        "HERC2_result": "blue",
        "HERC2_genotype": "GG",
    }


def test_merge_classifier_result_dict_without_result_key():
    results = {}
    _merge_classifier_result(results, "HERC2", {"genotype": "AG"})

    assert results == {"HERC2_genotype": "AG"}
