import typer


app = typer.Typer(
    name="molop",
    help="MolOP: Molecule OPerator CLI",
    add_completion=False,
    no_args_is_help=True,
)


@app.callback()
def main():
    """
    MolOP: Molecule OPerator CLI
    """
    pass


if __name__ == "__main__":
    app()
