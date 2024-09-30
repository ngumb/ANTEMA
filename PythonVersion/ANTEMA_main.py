""" Main script for Segmentation Analysis of HRTEM images
 Run this script to run the automated nanoparticle transmission electron
 micrograph analyzer (ANTEMA)


 Copyright (C) 2022  University of Duisburg-Essen
 
     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
 
     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


import glob
import matplotlib.pyplot as plt
from ncempy.io import dm
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import zscore
from skimage import morphology, measure
from tkinter import Tk, filedialog
import tensorflow as tf
from keras.utils import normalize
import os
import NPModel 





class ANTEMA:
    """ Main class for ANTEMA """

    def __init__(self,netname = 'NP-net-weights',Convexthresh = 0.95,PSep = True, saveas = 'xlsx', removecount = 0.5, minmarker = 0.5, SepMode = 'UECS'):
        """ Initialize with default parameters and add all inputs to self"""
        self.Convexthresh = Convexthresh
        self.PSep = PSep
        self.saveas = saveas
        self.removecount = removecount
        self.minmarker = minmarker
        self.SepMode = SepMode

    def batchmode(self):
        """Open file dialog to get the folder path of the folder containing all .dm3 files to be segmented and the paths to all the files"""
        root = Tk()
        root.attributes('-topmost',True)
        root.iconify()
        pathdir = filedialog.askdirectory(title = 'Select folder containing all .dm3 files to be segmented', parent = root)
        root.destroy()

        """Get all the files in the folder"""
        files = glob.glob(pathdir + '/*.dm3')
        """Ask How the files hould be saved, 0: one file for each image 1: all image data in one file 2: both"""
        self.saveas = int(input('How do you want to save the data?\n 0: one file for each image\n 1: all image data in one file\n 2: both'))
        """ Check if the input is 0,1 or 2 and if not repeat the question"""
        while self.saveas not in [0,1,2]:
            self.saveas = int(input('Input out of bounds. \nPlease type one of the following numbers:\n 0: one file for each image\n 1: all image data in one file\n 2: both'))
        for file in files:
            """Run the segmentation for each file"""
            self.run(file)
    
    def singlemode(self):
        """Open file dialog to get the path to the file to be segmented"""
        root = Tk()
        root.attributes('-topmost',True)
        root.iconify()
        path = filedialog.askopenfilename(title = 'Select dm3 file to be segmented', filetypes=[("dm3 files", "*.dm3")], parent = root)
        root.destroy()
        
        """Run the segmentation for the file"""
        self.run(path)
    
    def run(self, path):
        """Run the segmentation for the file"""
        
        """Load the dm3 file as a numpy array"""
        dm3 = dm.dmReader(path)
        image = dm3['data']
        pxSize = dm3['pixelSize'][0] 

        plt.figure()
        plt.imshow(image,cmap = 'gray')
        plt.show()
        """Apply a gauss filter and convert image to a uint8 rgb image"""
        image = gaussian_filter(image,2)
        '''normalize the image to values between 0 and 1 similarly to mat2gray in MATLAB'''
        image = self.mat2gray(image)
        image = (image*255).round().astype(np.uint8)

        '''apply zscore normalization with the values specified in the MATLAB code'''
        mean = 114.2725
        std =  58.9159
        image = (image - mean)/std
        image = np.stack((image,image,image),axis = 2)


        """Load the model"""
        model = NPModel.load_model([None,None,3], debug=True)

        '''Convert image to tensor and expand dimensions to fit the model input'''
        image = np.expand_dims(image, axis=0)
        image =  tf.convert_to_tensor(image, dtype=tf.float32)


        """Run the segmentation"""
        M = model.predict(image)
        
        M = np.squeeze(M)
        Mbackground = M[:,:,0]>0.5
        Mparticle = M[:,:,1]>0.5
        image = np.squeeze(image)
       
        exit()
        ### Code below is not finished and needs to be implemented
        """ Get particle properties"""
        props = self.particlePropertiesEval(Mparticle,pxSize,image )
        # Ab hier heavy testing!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        """Make image overlay"""
        overlay = self.overlay(image, M) ### probably already exists

        """ Add a scalebar to the image"""
        image = self.scalebar(image,pxSize)

        """Save the data""" #### not finished, needs to be implemented
        if self.saveas == 0:
            """Save each image in a seperate file"""
            pass
        elif self.saveas == 1:
            """Save all images in one file"""
            pass
        elif self.saveas == 2:
            """Save each image in a seperate file and save all images in one file"""
            pass
            
    def scalebar(self,image,pxSize):
        """Add a scalebar to the image"""
        ### not yet immplemented
        pass

    def particlePropertiesEval(self,M,pxSize,image):
        """Measure the particle properties"""
        """Make M into a integer array"""
        # Not yet implemented
        M = M.astype(int)
               
        """Remove all Areas with an equivalence diameter of lower than the set removecount diameter nm size"""
        removecount = int(np.pi/4*(self.removecount/pxSize)**2)
        M = morphology.remove_small_objects(M,removecount)

        """If separation routine is activated, separate particles with convexity lower than the threshold"""
        if self.PSep:
            """ Measure the particle properties
            Measure the boundingbox, perimeter, convex image and euler number of each particle"""
            label_image = measure.label(M)
            props = measure.regionprops(label_image, property = ['bbox','perimeter','image','image_convex','euler_number'])
            ''' Measure the convexity of each particle by dividing the perimeter of the particle by the perimeter of the convex image'''
            for prop in props:
                prop['convexity'] = prop['perimeter']/measure.perimeter(prop['image_convex'])

                if prop['convexity'] < self.Convexthresh:
                    """"Check if Euler Number is below -10 and if it is it should be removed from the evaluation and the image"""
                    if prop['euler_number'] < -10:
                        M[prop['image']] = 0 # CHECK IF THIS WORKS
                        """Give the user a warning message that the particle was removed from the image due to too many holes"""
                        Warning('Particle with too many holes removed from the image')
                    """Check which separation method to use"""
                    if self.SepMode == 'UECS':
                        Msplit = self.ParticleSeparationUECS(prop['image'],pxSize)
                        """Remove the old particle from the image and add the new particles as layers in the stack"""
                    elif self.SepMode == 'Watershed':
                        Msplit = self.ParticleSeparationWatershed(prop['image'],pxSize)
                        """Remove the old particle from the image and add the new particles as layers in the stack"""
        
        else:
            """Get the particle properties (area, equivalent diameter, perimeter, minimum and maximum feret diameter, boundinb box, centroid position and convexity )"""
            label_image = measure.label(M)
            props = measure.regionprops_table(label_image, property = ['area','perimeter','minor_axis_length','major_axis_length','bbox','centroid'])


        

    def ParticleSeparationUECS(self, M, pxSize):
        """Separate particles that are connected by a small bridge"""
        # Not yet implemented
        minmarker = int(np.pi/4*(self.minmarker/pxSize)**2)

        updated = True
        iter_count = 0
        while updated:
            updated = False
            label_image = measure.label(M)
            """Get convex area, area, perimeter, convex image and minor axis length of each particle"""
            props = measure.regionprops(label_image, property = ['area','perimeter','image','image_convex','minor_axis_length'])
            if np.mod(iter_count,2) == 1:
                se = morphology.disk(1)
            else:
                """ se2    = strel('arbitrary', ones(2,2));"""
                se = morphology.arbitary(2) # CHECK IF THIS WORKS
            
            for prop in props:


                iter_count += 1



        return Msplit
        pass

    def ParticleSeparationWatershed(self, M, pxSize):
        """Separate particles using the watershed algorithm"""
        ## Not yet implemented
        pass

    def mat2gray(self,image):
        """Normalize the image to values between 0 and 1 similarly to mat2gray in MATLAB"""
        image = (image - np.min(image))/(np.max(image)-np.min(image))
        return image

A =ANTEMA()
A.singlemode()

        









        

        




        




    



