from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='vpm',
      version='0.01',
      description='Vortex Panel Method',
      long_description="""Adapted for Python from Kuethe and Chow
      url='https://github.com/OwenShepherd/vortex-panel-method""",
      author='KCO',
      author_email='owen.shepherd@colorado.edu',
      license='MIT',
      packages=['vpm'],
      zip_safe=False,
      include_package_data=True)
