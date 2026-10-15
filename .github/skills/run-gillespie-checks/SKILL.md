---
name: run-gillespie-checks
description: Run focused or full tests for the Gillespie simulator and explain failures in the context of the current architecture.
disable-model-invocation: true
argument-hint: "Optional: provide a test path or change area"
---
Run the smallest relevant checks for the requested area, beginning with:

```bash
python -m pytest tests/gillespie -q
```

If a narrower path is supplied, use it first. Report passing checks, failures,
and whether each failure is caused by the current request or is pre-existing.

Use the user's request to choose a narrower test path when one is specified.
