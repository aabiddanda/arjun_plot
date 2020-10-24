import matplotlib.pyplot as plt

"""
  Various plotting functions
"""
def plot_yx(ax, **kwargs):
    """Plot the y=x line within a matplotlib axis
    """
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    # now plot both limits against each other
    ax.plot(lims, lims, **kwargs)

def debox(ax):
  """Remove the top and right spines of a plot
  """
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)


def label_multipanel(axs, labels=['A','B'], xoffset=-0.05, yoffset=1.14, **kwargs):
  """Labeling multiple axes with text labels
  """
  # Make sure that we have the nice labels here
  assert(len(axs) == len(labels))
  for i, lbl in enumerate(labels):
    axs[i].text(xoffset, yoffset, lbl, transform=axs[i].transAxes, **kwargs)
