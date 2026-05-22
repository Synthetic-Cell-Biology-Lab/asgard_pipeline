#!/usr/bin/env python3
"""Norns CLI for creating ASGARD pipeline configs."""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Optional

import typer
import yaml

app = typer.Typer(help="Norns: config creation CLI for ASGARD pipelines.")

# Paths are resolved relative to this script's location, not the caller's
# working directory. This is intentional: Norns is a monorepo-local tool and
# templates/processes live alongside it regardless of where it is invoked from.
REPO_ROOT = Path(__file__).resolve().parent.parent
TEMPLATES_DIR = REPO_ROOT / "templates" / "configs"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "processes"


def _template_name_to_path(template_name: str) -> Path:
    candidate = TEMPLATES_DIR / f"{template_name}.template.yaml"
    if not candidate.exists():
        # _list_templates() already returns a sorted list; no need for sorted() again
        available = ", ".join(_list_templates())
        raise typer.BadParameter(
            f"Unknown template '{template_name}'. Available templates: {available}"
        )
    return candidate


def _list_templates() -> list[str]:
    if not TEMPLATES_DIR.exists():
        return []
    return [
        path.name.replace(".template.yaml", "")
        for path in sorted(TEMPLATES_DIR.glob("*.template.yaml"))
    ]


@app.command("list-templates")
def list_templates() -> None:
    """List available config templates."""
    templates = _list_templates()
    if not templates:
        typer.echo("No templates found.")
        raise typer.Exit(code=1)
    typer.echo("Available templates:")
    for name in templates:
        typer.echo(f"- {name}")


@app.command()
def init(
    template: str = typer.Argument(..., help="Template name (without .template.yaml)."),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output path for generated config. Defaults to processes/<name>.yaml",
    ),
    name: Optional[str] = typer.Option(
        None,
        "--name",
        "-n",
        help="Base filename (without extension) when --output is not provided.",
    ),
    force: bool = typer.Option(False, "--force", help="Overwrite output if it exists."),
) -> None:
    """Create a new config file from a pipeline template."""
    # Warn when --name is redundant because --output already pins the path
    if output is not None and name is not None:
        typer.echo("Warning: --name is ignored when --output is provided.", err=True)

    template_path = _template_name_to_path(template)

    if output is None:
        stem = name or f"{template}_new"
        output = DEFAULT_OUTPUT_DIR / f"{stem}.yaml"

    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    if output.exists() and not force:
        raise typer.BadParameter(
            f"Output file already exists: {output}. Use --force to overwrite."
        )

    shutil.copyfile(template_path, output)
    typer.echo(f"Created config from template '{template}': {output}")


@app.command()
def validate(
    config: Path = typer.Argument(..., exists=True, dir_okay=False, readable=True)
) -> None:
    """Validate that a config file is parseable YAML and has a pipeline field."""
    # Catch malformed YAML before it produces an ugly unhandled traceback
    try:
        with config.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        typer.echo(f"Invalid YAML: {exc}")
        raise typer.Exit(code=1)

    if not isinstance(data, dict):
        typer.echo("Config is valid YAML but must be a mapping/object at top-level.")
        raise typer.Exit(code=1)

    # Treat pipeline: null and pipeline: "" as missing, not valid
    if not data.get("pipeline"):
        typer.echo("Config is missing required key: pipeline")
        raise typer.Exit(code=1)

    typer.echo(f"Config valid: {config} (pipeline={data['pipeline']})")


if __name__ == "__main__":
    app()