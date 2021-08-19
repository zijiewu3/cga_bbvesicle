from setuptools import setup

setup(
        name="cga_bbvesicle",
        install_requires="crease_ga",
        entry_points={"crease_ga.plugins":["bbvesicle=cga_bbvesicle.scatterer_generator:scatterer_generator"]},
        py_modules=["cga_bbvesicle"],
                        )
