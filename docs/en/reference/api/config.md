# molop.config

This page provides the API reference for the `molop.config` module.

The module exposes the process-wide `molopconfig` object. It controls progress
display, the maximum parallel-job limit, logging, and structure-recovery
options. Prefer explicit function arguments for per-call behavior; change the
global configuration only when the policy should apply to subsequent calls.

::: molop.config
