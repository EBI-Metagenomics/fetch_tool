import os
import sys


def pytest_generate_tests(metafunc):
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    os.environ["FETCH_TOOL_CONFIG"] = os.path.join(
        root, "config/fetchdata-config-default.json"
    )
    # Add fetchtool to the path
    sys.path.insert(0, root)
    os.environ["PATH"] += ":" + os.path.join(root, "fetchtool")
