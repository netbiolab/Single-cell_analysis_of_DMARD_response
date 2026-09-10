"""Minimal local replacement for training telemetry; never connects to a service."""
import argparse

config = None


class Settings:
    def __init__(self, **kwargs):
        pass


class Run:
    def log_artifact(self, *args, **kwargs):
        pass

    def finish(self):
        pass


class Artifact:
    def __init__(self, *args, **kwargs):
        pass

    def add_file(self, *args, **kwargs):
        pass


def init(config=None, **kwargs):
    globals()['config'] = argparse.Namespace(**config)
    return Run()


def log(metrics):
    # The original scGPT logger retains epoch losses in run.log.
    pass


def watch(*args, **kwargs):
    pass


def define_metric(*args, **kwargs):
    pass


def finish():
    pass
