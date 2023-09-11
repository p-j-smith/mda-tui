# MDAnalysis TUI
A textual user interface (TUI) for MDAnalysis trajectory transformations

https://github.com/p-j-smith/mda-tui/assets/29753790/b73c1a7e-ca97-41cd-b7a2-09899ecf1cb3

[//]: # (Badges)

| **Latest release** | [![Last release tag][badge_release]][url_latest_release] ![GitHub commits since latest release (by date) for a branch][badge_commits_since]  [![Documentation Status][badge_docs]][url_docs]|
| :----------------- | :------- |
| **Status**         | [![Tests][badge_test]][url_test] [![codecov][badge_codecov]][url_codecov] [![Linting][badge_linting]][url_linting] |
| **Community**      | [![License: GPL v2][badge_license]][url_license]  [![Powered by MDAnalysis][badge_mda]][url_mda]|

## Warning

MDanalysis TUI is alpha software and may change without warning.

[badge_test]: https://github.com/p-j-smith/mda-tui/actions/workflows/test_and_build.yaml/badge.svg
[badge_linting]: https://github.com/p-j-smith/mda-tui/actions/workflows/linting.yml/badge.svg
[badge_codecov]: https://codecov.io/gh/p-j-smith/mda-tui/branch/main/graph/badge.svg
[badge_commits_since]: https://img.shields.io/github/commits-since/p-j-smith/mda-tui/latest
[badge_docs]: https://github.com/p-j-smith/mda-tui/actions/workflows/docs.yml/badge.svg
[badge_license]: https://img.shields.io/badge/License-GPLv2-blue.svg
[badge_mda]: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
[badge_release]: https://img.shields.io/github/release-pre/p-j-smith/mda-tui.svg
[url_test]: https://github.com/p-j-smith/mda-tui/actions/workflows/test_and_build.yaml
[url_linting]: https://github.com/p-j-smith/mda-tui/actions/workflows/linting.yml
[url_codecov]: https://codecov.io/gh/p-j-smith/mda-tui/branch/main
[url_docs]: https://p-j-smith.github.io/mda-tui/
[url_latest_release]: https://github.com/p-j-smith/mda-tui/releases
[url_license]: https://www.gnu.org/licenses/gpl-2.0
[url_mda]: https://www.mdanalysis.org


## Installation

When installing MDAnalysis TUI, we highly recommend using virtual environments.
If possible, we strongly recommend that you use
[Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) as your package manager.

First ensure that you have [mamba](https://mamba.readthedocs.io/en/latest/index.html) installed.

Create a virtual environment and activate it:

```
mamba create --name mda-tui -c conda-forge python=3.10 pip
mamba activate mda-tui
```

Then install MDAnalysis TUI:

```
python -m pip install mda-tui
```

And when you are finished, you can exit the virtual environment with:

```
mamba deactivate
```

To use conda instead, replace all `mamba` commands with `conda`.


### Installation from source

To build MDAnalysis TUI from source, follow the above steps for creating a virtual environment, clone this repository:

```
git clone https://github.com/p-j-smith/mda-tui.git
```

Then inside your virtual environment type:

```
python -m pip install .
```

If you want to create a development environment, you can
install MDAnalysis TUI in editible mode along with
the dependencies required for running tests and and building docs with
the following command:

```
python -m pip install -e ".[dev]"
```

## Code of conduct

MDAnalysis TUI is bound by a [Code of Conduct](https://github.com/p-j-smith/mda-tui/blob/main/CODE_OF_CONDUCT.md).

## Copyright

The MDAnalysis TUI source code is hosted at https://github.com/p-j-smith/mda-tui
and is available under the GNU General Public License, version 2 (see the file [LICENSE](https://github.com/p-j-smith/mda-tui/blob/main/LICENSE)).

Copyright (c) 2023, Paul Smith

## Acknowledgements

Project based on the
[MDAnalysis Cookiecutter](https://github.com/MDAnalysis/cookiecutter-mda) version 0.1.
Please cite [MDAnalysis](https://github.com/MDAnalysis/mdanalysis#citation) when using MDAnalysis TUI in published work.
