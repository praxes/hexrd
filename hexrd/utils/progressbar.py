from __future__ import absolute_import

try:
    from progressbar import ProgressBar, Bar, ETA, Percentage, ReverseBar
except:
    # Dummy no-op progress bar to simplify code using ProgressBar
    class ProgressBar(object):
        def __init__(*args, **kwargs):
            pass

        def start(self):
            return self

        def finish(self):
            pass

        def update(self, x):
            pass

    Bar = ETA = Percentage = ReverseBar = ProgressBar
