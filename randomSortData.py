#Take in two lists (x-data and y-data for some points) and return
# the same lists, with a random rearrangement of which x-y pairings.

import random
import numpy as np

#This randomly shuffles a list of lists, keeping the pairings of list elements the same
def randomShuffleListOfLists(list_of_lists):
    zipped_list_of_lists = list(zip(*list_of_lists))
    random.shuffle(zipped_list_of_lists)
    new_list_of_lists = list(zip(*zipped_list_of_lists))
    new_list_of_lists = [list(new_arr) for new_arr in new_list_of_lists] 

    return new_list_of_lists


def randomSortList(list_to_randomize, list_to_randomize_errors = None, n_draws = -1, replace = 0):
    if n_draws < 0: n_draws = len(list_to_randomize)
    if not(replace):
        if list_to_randomize_errors is None: list_to_randomize_errors = [0 for elem in list_to_randomize]
        #print ('[list_to_randomize, list_to_randomize_errors] = ' + str([list_to_randomize, list_to_randomize_errors]))
        zip_list_with_errs = list(zip(list_to_randomize, list_to_randomize_errors))
        #print ('zip_list_with_errs  = ' + str(zip_list_with_errs ))
        random.shuffle(zip_list_with_errs)
        #print ('zip_list_with_errs  = ' + str(zip_list_with_errs ))
        new_list, new_list_errors = zip(*(zip_list_with_errs[0:n_draws]))

    else:
        new_list = [0 for elem in list_to_randomize]
        #if not list_to_randomize_errors is None: new_list_errs = [0 for elem in list_to_randomize]
        new_list_errors = [0 for elem in list_to_randomize]

        for i in range(n_draws):
            index = random.randint(0, len(list_to_randomize) - 1)
            #print 'x_index = ' + str(x_index)
            #print 'y_index = ' + str(y_index)
            new_list[i] = list_to_randomize[index]
            if not list_to_randomize_errors is None: new_list_errors[i] = list_to_randomize_errors[index]

    return [new_list, new_list_errors]

def randomSortListWithNans(list_to_randomize, list_to_randomize_errors = None, replace = 0):
    if list_to_randomize_errors is None: list_to_randomize_errors = [0 for elem in list_to_randomize]
    print ('Trimming nans...')
    trimmed_list = np.array(list_to_randomize)[np.logical_not(np.isnan(list_to_randomize))]
    flat_trimmed_errs = np.array(list_to_randomize_errors)[np.logical_not(np.isnan(list_to_randomize_errors))]

    print ('Sorting trimmed list...')
    sorted_trimmed_list, sorted_trimmed_list_errs = randomSortList(trimmed_list, flat_trimmed_errs, n_draws = -1, replace = replace)

    sorted_list = [0 for elem in list_to_randomize]
    sorted_errs = [0 for elem in list_to_randomize_errors]

    trim_index = 0
    print ('Rebuilding full list...')
    for i in range(len(sorted_list)):
        #if i % 100000 == 0: print ('i = ' + str(i))
        if np.isnan(list_to_randomize[i]):
            sorted_list[i] = np.nan
            sorted_errs[i] = np.nan
        else:
            sorted_list[i] = sorted_trimmed_list[trim_index]
            sorted_errs[i] = sorted_trimmed_list_errs[trim_index]
            trim_index = trim_index + 1

    return [ sorted_list, sorted_errs ]




def randomSortData(x_data, y_data, n_draws = -1, y_errs = None, replace = 0):

    if n_draws < 0: n_draws = len(x_data)

    new_x_data = []
    new_y_data = []
    if not y_errs is None: new_y_errs = []

    for i in range(n_draws):
        #print 'i = ' + str(i)
        x_index = random.randint(0, len(x_data) - 1)
        y_index = random.randint(0, len(x_data) - 1)
        #print 'x_index = ' + str(x_index)
        #print 'y_index = ' + str(y_index)
        new_x_data = new_x_data + [x_data[x_index]]
        new_y_data = new_y_data + [y_data[y_index]]
        if not y_errs is None: new_y_errs = new_y_errs + [y_errs[y_index]]

        if not replace:
            x_data = x_data[0:x_index] + x_data[x_index + 1:len(x_data)] if x_index < len(x_data) - 1 else x_data[0:len(x_data)-1]
            y_data = y_data[0:y_index] + y_data[y_index + 1:len(y_data)] if y_index < len(y_data) - 1 else y_data[0:len(y_data)-1]
            if not y_errs is None: y_errs = y_errs[0:y_index] + y_errs[y_index + 1:len(y_errs)] if y_index < len(y_errs) - 1 else y_errs[0:len(y_errs)-1]

    if not y_errs is None: return new_x_data, new_y_data, new_y_errs
    else: return new_x_data, new_y_data
