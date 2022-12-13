import setuptools
import glob
import HCGB.config.setup_module as setup_module

long_description_text = ""
with open("README.md", "r") as fh:
    long_description_text = fh.read()

setuptools.setup(
    name="RSP",
    version=setup_module.get_version("./VERSION"),

    scripts=glob.glob('main/*'),
    author="Jose F. Sanchez-Herrero, Mireia Marín Ginestar",

    author_email="jfbioinformatics@gmail.com",
    description="RNAseq pipeline",

    long_description_content_type="text/markdown",
    long_description=long_description_text,
    url="https://github.com/HCGB-IGTP/RNA_seq_pipe",
    packages=setuptools.find_packages(),
    license='MIT License',

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    install_requires=setup_module.get_require_modules("RSP/config/python/python_requirement_summary.txt"),
)
