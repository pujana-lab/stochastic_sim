---
name: Update project architecture
description: "Trace a requested project change and update the maintained architecture documentation without copying stale implementation notes."
argument-hint: "Describe the change or files whose architecture should be refreshed"
agent: "Architecture Maintainer"
---

Trace the requested change through source, tests, configuration adapters, and
user-facing entry points. Update
[`docs/ARCHITECTURE.md`](../../docs/ARCHITECTURE.md) only where the actual flow
or contract changed. Also update the documentation index or README when a
discoverability link is needed. Report contradictory historical documentation
without treating it as current behaviour.

Use the user's request as the change to trace.
