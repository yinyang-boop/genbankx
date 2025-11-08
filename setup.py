from setuptools import setup, find_packages

setup(
    name="genbankx",
    version="0.1.0",
    description="Lightweight GenBank parser with ordered CDS IntervalList and translation",
    author="BIO0001",
    license="MIT",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
)
