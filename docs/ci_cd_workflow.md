# CI/CD Workflow Overview

This document outlines the GitHub Actions workflows that power continuous integration and delivery for Pioneer.jl.

## Overview

The repository leverages several workflows to test, build, and release the project across platforms. The table below summarizes which events trigger each workflow.

| Workflow | Event Triggers | Purpose |
| --- | --- | --- |
| `tests.yml` | Push to any branch or tags matching `v*.*.*`; pull requests; manual dispatch | Run the test suite on Ubuntu with Julia 1.11 |
| `docs.yml` | Push to `main`; tags matching `v*`; pull requests; manual dispatch | Build and deploy documentation |
| `build_app_linux.yml` | Push to `main` or `develop`; tags matching `v*.*.*`; pull requests to `main` or `develop`; manual dispatch (optional `tag` input) | Build and package the Linux application |
| `build_app_macos.yml` | Push to `main` or `develop`; tags matching `v*.*.*`; pull requests to `main` or `develop`; manual dispatch (optional `tag` input) | Build and package the macOS application |
| `build_app_windows.yml` | Push to `main` or `develop`; tags matching `v*.*.*`; pull requests to `main` or `develop`; manual dispatch (optional `tag` input) | Build and package the Windows application |
| `CompatHelper.yml` | Scheduled daily (`0 0 * * *`); manual dispatch | Update package compatibility constraints |
| `TagBot.yml` | Tag pushes; creation of issue comments; manual dispatch | Tag new releases and interact with Julia's package registry |

## Additional Notes

- Most workflows support `workflow_dispatch` for manual invocation.
- Build workflows use concurrency controls to cancel in-progress runs when new commits arrive on the same reference.
- Tag-based releases (`v*.*.*`) trigger build and test workflows as well as documentation and tagging automation.

