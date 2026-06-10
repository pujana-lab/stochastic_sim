from pathlib import Path
import json as _json


def load_config_json(path: Path) -> dict:
    with open(path) as f:
        return _json.load(f)
