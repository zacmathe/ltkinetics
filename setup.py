from setuptools import setup, find_packages

setup(
    name='ltkinetics',
    version='0.3',
    packages=find_packages(exclude=['tests*', 'examples*']),
    license='GPL',
    description='simulate nitrogenase kinetics',
    long_description=open('README.md').read(),
    install_requires=['numpy', 'scipy'],
    url='https://github.com/zacmathe/ltkinetics',
    author='Zachary Mathe',
    author_email='zachary.mathe@cec.mpg.de',
)
