from setuptools import setup, find_packages

setup(
    name="probe-tester",
    version="1.0.0",
    description="Genomic primer/probe validation and reporting pipeline",
    author="Erick Samera",
    author_email="erick.samera@kpu.ca",
    packages=find_packages(),
    py_modules=["main"],
    install_requires=[
        "biopython",
        "rich",      # optional, but recommended
        "tabulate",  # optional, but recommended
    ],
    python_requires='>=3.8',
    entry_points={
        "console_scripts": [
            "probe-tester=main:main",
        ]
    },
    include_package_data=True,
)
