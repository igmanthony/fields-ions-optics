'''
Adapted in 2021 by Ian Anthony from MATLAB code written by Matthew Brantley
A companion script to https://igmanthony.com/fields-ions-optics/
This script does not implement some of the features found in the Simulation
Playground of the Simulating Fields, Ions, and Optics website and is
optimized for simplicity and readability, rather than performance. However,
the steps this script uses are almost identical, for static lenses, in the
generation and flight of ions, so results should be largely comparable.

For information on using or reading this script, please consult the
"readme.pdf" file or the page at:
https://igmanthony.com/fields-ions-optics/python_simulation_readme.html
'''

# We need to import some Python libraries to expand our capabilities!
# These first two libraries are included with Python itself
import os  # "Operating system" library for getting file paths
import math  # "Math" library for doing some mathematics functions

# These libraries are not always included, so you might need to download them
# Instructions for downloading and installing these libraries can be found in
# The "readme.pdf" file that accompanies this script
import cv2  # "Open CV" library for image loading and processing
# Note: Importing a library "as" something renames the library - so when you
# see "np" in the code, you know that we are working with the "numpy" library
import numpy as np  # "NumericalPython" library for doing matrix math
import matplotlib.pyplot as plt  # Plotting library for making plots


# "start()" is a function. Functions often have inputs and outputs. This
# "start" function does not have any inputs, so nothing is between the
# parentheses. "start" also does not have any outputs. Outputs are whatever
# comes after the "return" keyword, usually at the end of a function.
def start():
    '''
    The start of this script. "start" coordinates the simulation; it calls
    other functions and:
        - Assigns volages to the different-colored electrodes
        - Makes ions
        - Loads the png image of electrodes ("electrodes.png")
        - Generates the electric field from the png image of the electrodes
        - Flies the ions in the electric field
        - Plots the resulting ion flight paths
    '''

    # Let's load the 'electrodes.png' image and use it to make our electrode
    # and electric field maps! Because the computer does not know where our
    # image is, we need to orient ourselves and tell it the full file path to
    # the "electrodes.png" image which should be in the same directory as this
    # script
    image_name = 'electrodes.png'
    script_path = os.path.abspath(__file__)  # Get the full path of this script
    script_directory = os.path.dirname(script_path)  # Directory of this script
    # We join the directory to the image name to get the full path of the image
    image_path = os.path.join(script_directory, image_name)
    # Now we load the image as a numpy array of RGB colors
    rgb_image = cv2.cvtColor(cv2.imread(image_path), cv2.COLOR_BGR2RGB)
    # We convert the image from RGB colors to integer colors. This simplifies
    # the job of turning the electrode colors into the electrode voltages.
    # If you want to know how it works, scroll down to the line that starts
    # "def rgb_image_to_integer_image(rgb_image):"
    integer_image = rgb_image_to_integer_image(rgb_image)

    # ========================== MODIFY ME! ===================================
    # rgb_to_integer takes a red, a green, and a blue color value from 0 to 255
    red_electrode = rgb_to_integer(255, 0, 0)
    green_electrode = rgb_to_integer(0, 255, 0)
    blue_electrode = rgb_to_integer(0, 0, 255)
    # Let's define voltages for all the colored electrodes in the image. We
    # will put the colors and voltages into a "dictionary" so that we can look
    # up a color and find what voltage it is supposed to be.
    # Note: if we don't define an electrode's color, or if we give an electrode
    # a voltage of 0, it will be electrically invisible!

    # The curly-brackets --> {} define a dictionary, with keys (electrode
    # colors) and corresponding values (electrode voltages)
    electrodes = {
        red_electrode: 80,  # a red electrode at 80 V
        green_electrode: 100,  # a green electrode at 100 V
        blue_electrode: 1,  # a blue electrode at 1 V
    }  # consider defining more electrodes and adding them to the PNG image!
    # Don't forget: Every color in the PNG (except white) needs an electrode
    # color and an entry in the "electrodes" dictionary with a voltage

    # Now we can define the ions we want to simulate. We put the ions in a
    # list that is denoted with square-brackets --> [] . Inside this list,
    # each ion is a dictionary ( denoted with curly-brackets --> {} ).
    # Each ion needs a mass, a charge (z), a velocity, and a position.
    # The velocity and position will be lists with [x, y] coordinates.

    # We'll simulate five ions with identical parameters, except for their
    # x-velocities
    ions_to_simulate = [
        {'mass': 80, 'z': 1, 'position': [100, 195], 'velocity': [-250, 1]},
        {'mass': 80, 'z': 1, 'position': [100, 195], 'velocity': [-100, 1]},
        {'mass': 80, 'z': 1, 'position': [100, 195], 'velocity': [0.5, 1]},
        {'mass': 80, 'z': 1, 'position': [100, 195], 'velocity': [100, 1]},
        {'mass': 80, 'z': 1, 'position': [100, 195], 'velocity': [250, 1]},
    ]
    # =========================================================================

    # We turn the PNG image into a 'map' of the electrodes by replacing the
    # colors with voltages defined in the "electrodes" dictionary above
    electrode_map = make_electrode_map(integer_image, electrodes)
    # We generate electric fields out of the electrodes that we made. We will
    # use these fields to determine the electric gradient (and thus forces)
    # for the ions we will simulate
    # Note: this step sometimes takes a bit of time, so we'll leave ourselves
    # a message to let us know we've reached this point in the code
    print("Making the electric fields. Please wait...")
    electric_field_map = make_electric_field_map(electrode_map)

    # To visualize the ions flight paths, we should first set up a figure and
    # add to it as we simulate each ion. The left_plot will be the electrode
    # png image. The right_plot will be the electric field map. On both images,
    # we will plot the ions' flight paths.
    figure, (left_image, right_image) = plt.subplots(1, 2)

    # We can loop through multiple ions to see how varying a parameter changes
    # the simulation. Note: This is a "for" loop: it can be read as "for each
    # ion in the list of ions to simulate, do:"
    for ion in ions_to_simulate:
        # Now we're in the loop (indented by 4 spaces). So let's fly the ion!
        # Read down to the "def fly_ion" function to see what it does!
        print("Flying an ion. Please wait...")  # this can take some time
        xs, ys = fly_ion(ion, electrode_map, electric_field_map)
        # now let's plot the resulting flight path on both images!
        left_image.plot(xs, ys)
        right_image.plot(xs, ys)

    # Now we've completed the loop (notice we de-indented by 4 spaces)
    # These next few lines of code are not needed for understanding the
    # simulation. They are for making the plots prettier.
    left_image.axis('off')  # removing the axis makes the plot a bit cleaner
    left_image.set_title('Ion paths', fontsize=18)  # add a title
    left_image.imshow(rgb_image)  # show the RGB image (png) we loaded with cv2
    right_image.axis('off')  # removing the axis makes the plot a bit cleaner
    right_image.set_title('Electric field', fontsize=18)  # add a title here
    color = right_image.imshow(electric_field_map, cmap='turbo')  # nice colors
    colorbar = figure.colorbar(color, fraction=0.046, pad=0.04)
    colorbar.set_label('Volts', rotation=90, fontsize=18)
    plt.show()  # now let's show the plot and print a completion message!
    print("\nSimulation completed successfully!")


def rgb_image_to_integer_image(rgb_image):
    '''transforms the r,g,b values of an image (3D numpy array) to integers'''
    rgb_image = np.array(
        rgb_image, dtype=int)  # change the datatype to integer
    # we extract each color component from the image
    r, g, b = rgb_image[:, :, 0], rgb_image[:, :, 1], rgb_image[:, :, 2]
    # then call the rgb_to_integer function for the arrays of colors
    integer = rgb_to_integer(r, g, b)
    return integer


def rgb_to_integer(red, green, blue):
    '''converts red, green, and blue numbers (each 0 to 255) to an integer'''
    # "<<" is a bit-shift operator. It moves the binary 0's and 1's to the
    # left. We use it to turn 3 small numbers into one larger number.
    integer = (((red << 8) + green) << 8) + blue
    return integer


def make_electrode_map(electrode_image, electrodes):
    '''Converts the image of colors into a map of voltage values
    where each pixel of color is replaced with the corresponding voltage.'''
    # we start by making a map of all zeros that we will add to
    electrode_map = np.zeros_like(electrode_image)
    # then we loop, using a "for" loop through all the colors and voltages we
    # defined in the "electrodes" dictionary in the "start" function
    # Read as "for each color and voltage in the "electrodes" dictionary do:"
    for color, voltage in electrodes.items():
        # Make the pixels with colors equal to "color" equal to the "voltage"
        electrode_map[electrode_image == color] = voltage
    return electrode_map


def make_electric_field_map(electrode_map):
    '''Generates an electric field for all electrodes by performing the finite
    difference method (FDM) to solve the Laplace equations'''
    # the electric field map (e_field_map) starts off as zeros - we'll add to
    # it to create the final map and return this variable at the end
    e_field_map = np.zeros_like(electrode_map, dtype=np.float64)
    # Read this as: "For each unique volts value in the electrode map do:"
    for volts in np.unique(electrode_map):  # quick way to iterate electrodes
        if volts == 0:  # the "white" color in the png is a "0 volt" electrode
            continue   # so we skip this color/volt value
        # We should make an image of just the electrode with the volts we are
        # wanting to look at as well as all the "other" electrodes
        this_mask = electrode_map == volts  # we mask off the current electrode
        other_mask = electrode_map != 0  # we also mask all electrodes
        other_mask[this_mask] = False   # and remove just the current electrode

        # We make a passive electric field of 0.5 volts (this is just a
        # low number that is not 0, it could be 0.3 volts and work the same).
        e_field = np.ones_like(electrode_map) / 2  # 1 / 2 = 0.5 volts
        # We apply the mask to the e_field to set it to our "volts" value
        e_field[this_mask] = volts
        # We apply the mas to all other e_fields and set them to 0 volts
        e_field[other_mask] = 0
        # We will downsample the image and the upscale it to make this a bit
        # faster (working on small images is faster).
        for division in [32, 16, 8, 4, 2, 1]:
            # perform an FDM step to propagate the electric field
            e_field = fdm_step(e_field, volts, this_mask, other_mask, division)
            # upscale to the original size with interpolation
            e_field = cv2.resize(
                e_field, electrode_map.shape, interpolation=cv2.INTER_LINEAR
            )
        # ========================= NOW YOU KNOW ==============================
        # Having a single electric field map is convenient, but assumes a
        # static field (no electrodes ever change). If we wanted a dynamic
        # field (for instance, a quadrupole), we would need to have some
        # function to calculate the electric field over and over at each time
        # point
        # Note: Be careful when using a dynamic field that the "phase" of the
        # electric fields are randomized (if they would be normally) - solutions
        # and simulations can be found for a set of deterministic changes to an
        # electric topology that would be impossible if the phases of the
        # electrodes were randomly distributed (as often happens in the real
        # world)
        # =====================================================================

        # now that we've calculated the electric field for an electrode, we add
        # that field to our total electric field map and move on to the next
        # electrode.
        e_field_map += e_field
    return e_field_map  # FYI: the e_field map is in volts per meter


def fdm_step(e_start, volts, this_mask, other_mask, div, conv=0.01):
    '''Calculate a single "step" using the finite difference method (FDM) to
    solve the laplace differential equation for the electric field gradient
    maps'''
    # duplicating the starting electric field is needed to prevent altering
    # the original data
    e_prev = np.copy(e_start)
    # some initial error is created by updating laplace boundaries
    e_now = update_laplace_boundaries(e_start, this_mask, other_mask, volts)
    # downsampling then upsampling allows for efficient field propagation
    # we do this for both "steps" of the electric field (e_now and e_prev) as
    # well as for the two masks of this electrode and the other electrodes
    # these four lines are downsampling the masks and fields
    e_now = e_now[::div, ::div]
    e_prev = e_prev[::div, ::div]
    this_mask = this_mask[::div, ::div]
    other_mask = other_mask[::div, ::div]
    # we need to know what the final row and column counts are to generate
    # a new array with a slightly different shape
    row, col = e_now.shape
    # calculate the error between the non-updated Laplace and updated Laplace
    # fields, we will loop and each loop should decrease this error we will
    # stop looping when the error is low enough (less than our "convergence"
    # or "conv") value defined in the "fdm_step" function signature above
    error = np.max(np.abs(e_now - e_prev))

    # Instead of a "for" loop, we use a "while" loop to loop until a condition
    # is met. In this case the condition is that the "error" value is less than
    # our convergence "conv" that we define in the function signature above.
    while error > conv: # we just loop until the error is acceptably low
        # We'll make a new (slightly larger) array of zeros and then paste
        # four offset copies of the e_now array into our new array
        e_new = np.zeros((row + 2, col + 2))
        e_new[1:-1, 0:-2] += e_now
        e_new[0:-2, 1:-1] += e_now
        e_new[2:, 1:-1] += e_now
        e_new[1:-1, 2:] += e_now
        # Now we'll divide by 4, so that we have "averaged" or blurred the
        # values of the cells of our array together.
        # ===================== NOW YOU KNOW ==========================
        # This is the finite difference method, learn more about this
        # method of generating electric fields from papers such as:
        # Dhumal, M. L.; Kiwne, S. B.; IJSAM, 9 (1), 2014 pp.11-13
        # =============================================================
        e_new[1:-1, 1:-1] /= 4

        # As we have eroded the boundaries of the electrodes, we need to
        # reinforce them and remove abberations at the edges of our simulation
        e_now = update_laplace_boundaries(
            e_new[1:-1, 1:-1], this_mask, other_mask, volts
        )
        # recalculate error, if it's too large the while loop will repeat
        error = np.max(np.abs(e_now - e_prev))
        e_prev = e_now.copy()  # update the "previous" value for next loop
    return e_now


def update_laplace_boundaries(e_field, this_mask, other_mask, volts):
    '''Ensures that the geometry is not eroded after the averaging and also
    correct external boundaries to appear infinite - this also means that
    any pixels around the edges of the .png image will dissappear!'''
    # correct external boundaries by copying the pixels next to the edges out
    # to the edges
    e_field[0, :] = e_field[1, :]
    e_field[-1, :] = e_field[-2, :]
    e_field[:, 0] = e_field[:, 1]
    e_field[:, -1] = e_field[:, -2]
    # reapply masks to prevent erosion of the electrodes (otherwise, it would)
    # appear as if the electrodes were turned off between steps
    e_field[this_mask] = volts
    e_field[other_mask] = 0
    return e_field


def fly_ion(ion, e_map, e_field_map):
    '''
    Calculates each time step (1 ns) of the ion simulation for up to 400,000
    time steps (400 microseconds total). The function will end if the ion
    collides or "splats" with an electrode or the simulation boundary. The
    function will return a set of x and y coordinates of the ion for each time
    step (as well as the time in microseconds). These data allow for plotting
    the ion flight trajectory, as well as performing other calculations (such
    as finding statistics about the ion's velocity). If the ion splats, the
    list of x, y points will be shorter than the maximum number of time steps
    and the last set of x, y points will be the splat location.

    Note: Mathematically, what we are doing in flying ions is similar to 
    moving a ball down a frictionless hill, if the topographical "map" of the
    "hill" is the same as our "e_field_map" or height-map of voltages.
    '''

    # Note: We use the "x,y" position and velocities for both pixel-values and
    # for the kinematics equations below. This means that 1 pixel will be
    # equivalent to 1 meter - a distance far too large for our simulation size.
    # (e.g. 200 x 200 pixels would be 200 x 200 meters)
    # So we will scale the delta time (dt) up 1,000-fold. Scaling time up
    # this much will make it appear as though the simulation has a pixel size:
    # 1 pixel = 1 mm
    time_scale = 1000
    # each step should be 1 ns (1e-9 s).
    dt = 1e-9 * time_scale
    # We will simulate for 400 microseconds (400_000) nanoseconds
    max_time = 400_000

    # We store some values locally and convert units to SI-standard
    position = np.array(ion['position'], dtype=float)
    velocity = np.array(ion['velocity'], dtype=float)
    mass = ion['mass'] * 1.660538921e-27  # convert mass to kg
    charge = ion['z'] * 1.602e-19  # convert charge to Coulombs
    x_boundary, y_boundary = e_map.shape
    x, y = position

    # the data we will save from this function are the x, and y values
    # of the ion as it moves. We'll begin these lists here and add to them
    # as we calculate each step of the ion's journey
    x_log, y_log = [x], [y]

    # Each step of this loop moves ahead by delta time (dt).
    # Each dt:
    #    the electric field is found for the ion's location
    #    the electric field is used to calculate the force on the ion
    #    the force on the ion is used to calculate the acceleration of the ion
    #    the acceleration of the ion is used to calculate the ion's velocity
    #    the velocity of the ion is used to calculate the ion's new position
    #    the new position of the ion is stored for the next timepoint
    #    the ion's new position is checked to see if the ion collides (splats)
    #        with anything (a wall or an electrode)
    #    if it has, then the function ends by returning the ion's times, x, and
    #        y values, otherwise the loop repeats 1 dt later! until 400_000 dt
    for step in range(1, max_time + 1):
        # first we get the electric vector at the location of our ion
        electric_vector = get_electric_vector(e_field_map, x, y)
        # let's calculate some simple kinematics parameters to find the new
        # ion position
        # ============================ NOW YOU KNOW ===========================
        # We are using numpy arrays of both x and y coordinates. Simple
        # mathematical expressions such as division and addition are propagated
        # to each element (in this case, the x and y values) of the numpy
        # array.
        # If we were not using numpy, we would need to double the number of
        # steps to mathematically calculate all parameters in both the x and y
        # directions (e.g.,
        #   force_x = charge * electric_vector_x
        #   force_y = charge * electric_vector_y
        # )
        # Adding a third dimension is somewhat easy, as it just means extending
        # these calculations to "z"s as well as xs and ys.
        # =====================================================================
        force = charge * electric_vector
        acceleration = force / mass
        velocity += acceleration * dt
        position += velocity * dt
        x, y = position  # update for logging and the next loop
        x_log.append(x)  # store the new x and y locations as well as the time!
        y_log.append(y)

        # check to see if the ion has splatted and the simulation should end
        if ion_has_splatted(int(y), int(x), e_map, x_boundary, y_boundary):
            # if the ion has splatted, then let's print ourselves a message
            # letting us know when and where it splatted
            print(f'SPLAT! {dt * step * 1_000:.3f} us; position: {position}')
            break
    else:
        # We don't want the simulation to go on forever, so we have to end
        # at the "maximum time". Let's let ourselves know the ion didn't splat.
        print('Maximum time reached:', dt * step * 1000, 'us', step, 'steps')
    return x_log, y_log


def get_electric_vector(e_field, y, x):
    '''finds the electric vector at the ion's location'''
    # Ugly math code
    ylow, xlow = int(y), int(x)
    yhigh, xhigh = min(math.ceil(y), 199), min(math.ceil(x), 199)
    yround, xround = min(int(round(y)), 199), min(int(round(x)), 199)

    # here we look up the voltage at 4 points and calculate the differences
    # between them to generate a "slope" in the x and y directions
    voltage_slope_y = (e_field[xhigh, yround] - e_field[xlow, yround])
    voltage_slope_x = (e_field[xround, yhigh] - e_field[xround, ylow])
    # invert the vector - as the ion is travelling opposite to the "hill" -
    # otherwise the ion would move towards the like-charged electrodes
    return -np.array([voltage_slope_x, voltage_slope_y])


def ion_has_splatted(ion_x, ion_y, e_map, x_boundary, y_boundary):
    '''checks whether the ion has hit any walls or an electrode'''
    # this is a chained "or" boolean conditional, so if anything is true, then
    # it is all true
    # Read this as "the ion has splatted if it hit a wall or an electrode"
    ion_splatted = (
        ion_x < 1                       # the ion hit the left-wall
        or ion_y < 1                    # the ion hit the top-wall
        or ion_x > x_boundary - 1       # the ion hit the right-wall
        or ion_y > y_boundary - 1       # the ion hit the bottom-wall
        or e_map[ion_x, ion_y] != 0     # the ion hit an electrode
    )
    return ion_splatted  # this "ion_splatted" value is either "True" or "False"


# This is where the program actually "starts" after the imports and "making"
# the functions - this calls the "start" function, back at the top of the
# program
start()
