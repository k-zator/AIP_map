import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="AIP_interaction_map", # Replace with your own username
    version="2.0.0",
    author="K J Zator; M C Storer",
    author_email="kz265@cam.ac.uk",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    package_data={'AIP_interaction_map': ['data/*', "mdtraj/*"], },
    entry_points={ 'console_scripts': ['Package = AIP_interaction_map.__main__:main' ]},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
