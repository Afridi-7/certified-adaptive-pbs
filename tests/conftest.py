"""Shared settings for the validity test suite (run: python -m pytest -q tests)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def pytest_sessionfinish(session, exitstatus):
    import json, platform
    import numpy, scipy
    from helpers import RESULTS
    out = os.path.join(os.path.dirname(__file__), "..", "results", "validity_tests_v1")
    os.makedirs(out, exist_ok=True)
    json.dump(dict(exitstatus=int(exitstatus), alpha=0.05, python=platform.python_version(), numpy=numpy.__version__,
                   scipy=scipy.__version__, results=RESULTS), open(os.path.join(out, "report.json"), "w"), indent=1)
