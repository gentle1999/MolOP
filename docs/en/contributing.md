# Contributing

Thank you for your interest in MolOP! We welcome contributions in all forms.

- **Reporting Bugs**: How to submit an effective Issue.
- **Feature Suggestions**: Share your ideas for improving MolOP.
- **Code Contribution Workflow**: Detailed steps from Fork to Pull Request.
- **Development Environment Setup**: How to configure your local development environment.
- **Code Style Guidelines**: Follow Ruff and type checking requirements.
- **Documentation Style Guide**: Use the unified four-layer documentation standard (`style_guide.md`).

## Documentation Quality Assurance

To ensure high-quality documentation, we enforce the following policies:

- **CI Verification**: Every Pull Request triggers a documentation build using `mkdocs build --strict`. This ensures no broken internal links and valid configuration.
- **Translation Policy**:
  - We use `TODO(translate):` as a placeholder for content that has not yet been translated between English and Chinese.
  - Placeholders are allowed in the `main` branch and Pull Requests to enable staged parity.
  - **Release Blockers**: Placeholders are strictly forbidden in release tags (`v*`). The CI will fail if any are detected during a release build.
- **Notebooks**: Jupyter notebooks in the documentation are optional and are not executed by the CI. Please ensure they are pre-executed if you want the output to be visible.
