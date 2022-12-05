import click

from .utils import Config
from . import processing
from . import explorer
from . import figures


@click.group(
    help='Monitor the pointing offset SPICE',
    chain=True,
    )
@click.option(
    '-c', '--conf-file',
    help='Configuration file',
    required=True,
    type=click.Path(exists=True),
    )
@click.pass_context
def cli(ctx, conf_file):
    ctx.obj = Config(conf_file)


@cli.command(help='Computes the offset between SPICE and FSI')
@click.pass_obj
def process(conf):
    processing.process(conf)


@cli.command(help='Generates an interactive result explorer')
@click.pass_obj
def gen_explorer(conf):
    explorer.gen_explorer(conf)


@cli.command(help='Generates the figures')
@click.pass_obj
def gen_figures(conf):
    figures.gen_figures(conf)
