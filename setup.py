from setuptools import find_packages, setup

setup(
    name="GravNN",
    packages=find_packages(),
    version="0.8.0",
    description="GravNN Package",
    author="John M",
    license="MIT",
    setup_requires=[
        "numpy == 1.23.0",
        "matplotlib",
        "numba",
        "pandas",
        'tensorflow>=2.6; sys_platform != "darwin"',
        'tensorflow-macos; sys_platform == "darwin"',
        "scikit-learn",
        "plotly",
        "dash",
        "trimesh == 3.9.31",
        "sigfig",
        "tqdm",
        "colorama",
        "tensorflow_model_optimization",
        # 'kaleido',
        "rtree",
        "pooch",
        "seaborn",
        "spicepy",
        "sphinx",
        "sphinx_rtd_theme",
        "sphinx-gallery",
        "pre-commit",
    ],
    install_requires=[
        "numpy == 1.23.0",
        "matplotlib",
        "numba",
        "pandas",
        'tensorflow>=2.6; sys_platform != "darwin"',
        'tensorflow-macos; sys_platform == "darwin"',
        "scikit-learn",
        "plotly",
        "dash",
        "trimesh == 3.9.31",
        "sigfig",
        "tqdm",
        "colorama",
        "tensorflow_model_optimization",
        # 'kaleido',
        "rtree",
        "pooch",
        "seaborn",
        "spicepy",
        "sphinx",
        "sphinx_rtd_theme",
        "sphinx-gallery",
        "pre-commit",
    ],
    tests_require=["pytest==4.4.1"],
    test_suite="tests",
)
