sint (Source INTensities estimation)
=================================================

 |Python Version|

|Tests| |Codecov|

|pre-commit| |Black|

.. |PyPI| image:: https://img.shields.io/pypi/v/funnel.svg
   :target: https://pypi.org/project/funnel/
   :alt: PyPI
.. |Python Version| image:: https://img.shields.io/pypi/pyversions/funnel
   :target: https://pypi.org/project/funnel
   :alt: Python Version
.. |Tests| image:: https://github.com/avivajpeyi/funnel/workflows/Tests/badge.svg
   :target: https://github.com/avivajpeyi/funnel/actions?workflow=Tests
   :alt: Tests
.. |Codecov| image:: https://codecov.io/gh/avivajpeyi/funnel/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/avivajpeyi/funnel
   :alt: Codecov
.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit
.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Black

Corresponding Paper
------------
https://www.overleaf.com/read/hytzvjcqxxqv#bfbf31

Contributing
------------

How to set up your development environment
------------------------------------------

You need Python 3.6+ and the following tools:

- Poetry_
- Nox_

Install the package with development requirements:

.. code:: console

   $ poetry install

You can now run an interactive Python session,
or the command-line interface:

.. code:: console

   $ poetry run python
   $ poetry run funnel

.. _Poetry: https://python-poetry.org/
.. _Nox: https://nox.thea.codes/


How to test the project
-----------------------

Run the full test suite:

.. code:: console

   $ nox

List the available Nox sessions:

.. code:: console

   $ nox --list-sessions

You can also run a specific Nox session.
For example, invoke the unit test suite like this:

.. code:: console

   $ nox --session=tests

Unit tests are located in the ``tests`` directory,
and are written using the pytest_ testing framework.

.. _pytest: https://pytest.readthedocs.io/

How to lint the project
-----------------------

You can ensure that your changes adhere to the code style by reformatting with Black_:

.. code:: console

   $ nox --session=black

It is recommended to open an issue before starting work on anything.
This will allow a chance to talk it over with the owners and validate your approach.

.. _pull request: https://github.com/avivajpeyi/funnel/pulls
.. _Black: https://black.readthedocs.io/
.. github-only
