from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
      name='aerotools',
      version='0.01',
      description='Aerospace Engineering Toolbox',
      long_description= readme(),
      author='Owen Shepherd',
      author_email='owen.shepherd@colorado.edu',
      license='MIT',
      packages=setuptools.find_packages(),
      zip_safe=False,
      include_package_data=True
      install_requires=[
          'numpy'
      ]
)
