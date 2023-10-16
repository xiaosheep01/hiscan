import setuptools

with open("README.md", "r",encoding="utf-8") as fh:
  long_description = fh.read()

setuptools.setup(
  name="hiscan",
  version="1.0",
  author="Yang Xiao",
  author_email="fredrik1999@163.com",
  description="Scanning histone mimics in secquences.",
  keywords="histone, mimic, viral protein",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/xiaosheep01/hiscan",
  packages=setuptools.find_packages(),
  install_requires=["colorama", "pandas", "numpy"],

  classifiers=[
        "Development Status :: 4 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Statistic",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
        "Operating System :: OS Independent",
    ],
  entry_points={
             'console_scripts': [
                 'hiscan = hiscan.main:starts',
             ],
    }
)
