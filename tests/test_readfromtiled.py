"""Unit tests for readfromtiled URI construction and result conversion.

No live Tiled server is needed — only the URI strings and pure conversion
logic are tested here.
"""

import pytest

pytest.importorskip("requests")

from matilda.readfromtiled import _build_search_uri, convert_results


def test_uri_basic_plan_filter():
    uri = _build_search_uri(10, plan_name="Flyscan")
    assert "page[limit]=10" in uri
    assert "filter[eq][condition][key]=plan_name" in uri
    assert 'filter[eq][condition][value]="Flyscan"' in uri
    assert "sort=-time" in uri
    assert "select_metadata={" in uri
    assert "filter[time_range]" not in uri
    # the exit_status FILTER is off by default (it does appear in
    # select_metadata, which is fine)
    assert "filter[contains][condition][exit_status]" not in uri


def test_uri_time_range():
    uri = _build_search_uri(5, plan_name="SAXS", lastNdays=2)
    assert "filter[time_range][condition][since]=" in uri
    assert "filter[time_range][condition][until]=" in uri
    assert "filter[time_range][condition][timezone]=US/Central" in uri


def test_uri_title_regex_and_exit_status():
    uri = _build_search_uri(3, plan_name="WAXS", title_regex="(?i)blank",
                            require_exit_status=True)
    assert "filter[contains][condition][exit_status]" in uri
    assert "filter[regex][condition][key]=title" in uri
    assert "filter[regex][condition][pattern]=(?i)blank" in uri


def test_uri_single_condition_per_filter_type():
    """Tiled ignores repeated filter keys — the builder must never emit two
    conditions of the same filter type."""
    uri = _build_search_uri(10, plan_name="Flyscan", title_regex="(?i)blank",
                            lastNdays=1, require_exit_status=True)
    assert uri.count("filter[eq][condition][key]") == 1
    assert uri.count("filter[regex][condition][key]") == 1


def test_uri_custom_select_metadata():
    uri = _build_search_uri(1, plan_name="tune_ar", select_md="uid:start.uid")
    assert uri.endswith("select_metadata={uid:start.uid}")


def _entry(uid, plan_name="Flyscan", exit_status="success",
           hdf5_file="f.h5", hdf5_path="/data"):
    return {
        "id": uid,
        "attributes": {"metadata": {"selected": {
            "plan_name": plan_name,
            "exit_status": exit_status,
            "hdf5_file": hdf5_file,
            "hdf5_path": hdf5_path,
        }}},
    }


def test_convert_results_uses_inline_exit_status():
    """exit_status from select_metadata avoids per-uid requests (N+1 fix);
    aborted runs and runs without an hdf5 file are skipped."""
    r = {"data": [
        _entry("u1"),
        _entry("u2", exit_status="abort"),
        _entry("u3", hdf5_file=None),
    ]}
    out = convert_results(r)
    assert out == [["/data", "f.h5"]]
