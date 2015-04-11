from setuptools import setup

setup(name='rebound',
    version='0.1',
    description='The N-body integrator REBOUND',
    url='http://github.com/hannorein/rebound',
    author='Hanno Rein',
    author_email='hanno@hanno-rein.de',
    license='GPL',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],
    keywords='astronomy astrophysics nbody integrator symplectic wisdom-holman',
    packages=['rebound'],
    install_requires=[],
    zip_safe=False)
