# Plasma-Spectroscopy-Gradient-Descent
Find out what's in your inductively coupled plasma using some calculus, linear algebra, and a dash of creativity.<p>
Here's a [paper](https://drive.google.com/file/d/1FzJek8TylvBXz0bf5_gJxKRb-Vu4l2ON/view?usp=sharing) I wrote on this project.</p>
## About This Project
This project was created to add another layer of quantification to the characterization of an oxygen plasma in an ultra-high-vacuum environment. Characterization of the plamsa is important because there are chemical reactions which lead to certain byproducts (ions) that can be extracted from the plasma for applications in nanofabrication, and by optimizing the conditions inside the plasma, we can maximize the amount of the extracted ions of interest. One way to probe these characteristics is to observe the light which is emitted due to excitational collisions between the gases and electrons inside the chamber. By identifying the wavelength of the emissive light, we can begin to learn about what types of ions, atoms, and molecules are present. However, this tells us little about the relative composition among the various species of gases present in the chamber, so there is motivation to find a more quantitative way of characterizing the emitted light.
<IMG SRC="spectrum3-9.gif">
### Start Here
The notebook will walk you through the experimental process, basic background, and strategies used.  
- [Jupyter Notebook](Spectroscopy%Gradient%Descent.ipynb)
  
After seeing how the process works, feel free to look through the source code or try your own characterization strategies using the data. If you have any ideas on how to improve on the groundwork I lay out here, I'm open to suggestions.
#### Most Useful Source Code:
- [Gradient Descent Algos](SpectrumFitting.py)  
- [Visualization Tools](GradVisualization.py)  
#### Supporting Code Base:
- [Reference Spectrum Handling/Processing](LineID.py)
- [Spectrum Data Processing](SpectralDataProcessing)
#### Data Files
- [Database Reference Lines](OFeNHlines.txt)
- [Spectrum Test Data](test_data.csv)
  
  
There's a lot more functionality in the code base than what is actually shown here, but the data used for that part of the project unfortunately can't be shared publicly.
