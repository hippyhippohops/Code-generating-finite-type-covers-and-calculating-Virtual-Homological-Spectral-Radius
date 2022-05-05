import itertools
from itertools import permutations
import networkx
import numpy
import numpy as np
from sympy.combinatorics.free_groups import free_group, vfree_group, xfree_group
import random
import networkx as nx
import matplotlib.pyplot as plt
from sympy import *
import pandas as pd
from tabulate import tabulate
from treelib import Node, Tree
from numpy import linalg as LA
import math


def permutation_list(size_of_list):
   list_n = []
   for i in range(1, size_of_list+1):
       list_n.append(i)
   s_n = permutations(list_n)
   return s_n

#Function to build set of permutation matrices given set of permutations
def build_permutation_matrices(s_n, size_of_list):
   permutation_matrices = []
   for i in s_n:
       matrix = np.zeros((size_of_list, size_of_list))
       for j in range(0, size_of_list):
           m = i[j]
           matrix[j][m-1] = 1
       permutation_matrices.append(matrix)
   return permutation_matrices

#function to generate free_group
def free_group_generator(generators_for_free_group):
   F = vfree_group(generators_for_free_group)
   return F

#Construct the representation from F to Permutation matrices
def representation_function(generators_for_free_group, matrix_1, matrix_2):
   n = len(generators_for_free_group)
   selected_matrices = [matrix_1, matrix_2]
   y = dict()
   j = 0
   for i in generators_for_free_group:
       if i != ",":
           y[i] = selected_matrices[j]
           j = j+1
   return y

#First, we construct the vertices of the covering space
def covering_space_vertices(size_of_list):
   n = size_of_list
   vertex_set = []
   for i in range(n): #Make sure indexing is okay
       a = np.eye(1, n, i)
       vertex_set.append(a)
   return vertex_set

#Function to construct the covering space
#Need to make sure matrices act on the right. Since matrices act on the right, ensure that row vector x matrice multiplication is correct
def generating_covering_space(generators_of_free_group, covering_vertex_set, rho):
   covering_space_edge_set = []
   for i in covering_vertex_set:
       for j in generators_of_free_group:
           if j != ",":
               a = np.array(i)
               b = rho.get(j)
               covering_space_edge_set.append([a,a.dot(b)])
   return covering_space_edge_set

def drawing_covering_space(covering_vertex_set, covering_space_edges, rho, generators_for_free_group, size_of_list):
   graph = nx.MultiDiGraph()

   #I will create this dictionary for the loops to become words later on
   f = dict() #dict_from_loops_to_words
   set_to_create_dict = []

   #Create vertices label
   y = dict()
   n = len(list(covering_vertex_set))
   i=0
   for j in range(1, n+1):
       y[j] = covering_vertex_set[i]
       i = i+1

   #Add nodes to graph
   for i in range(1,n+1):
       graph.add_node(i)
   #print(graph)

   #Add edges
   nodes_list = []
   for i in covering_vertex_set:
       for j in generators_for_free_group:
           if j != ",":
               a = np.array(i) #starting vertex of an edge
               b = rho.get(j) # b is Matrix for matrix mutliptication
               x = a.dot(b) #ending vertex of an edge
               for k in range(1,size_of_list+1):
                   if np.array_equal(y.get(k), a):
                       starting_vertex = k
                   if np.array_equal(y.get(k), x):
                       ending_vertex = k
               edge = (starting_vertex, ending_vertex)
               if edge in f.keys():
                   numpy.append(f[edge] , b)
               else:
                   f[edge] = b
               set_to_create_dict.append((edge,b))
               nodes_list.append(edge)
   graph.add_edges_from(nodes_list)
   return graph, nodes_list, f

#This function provides us the edges of the spanning tree

def create_spanning_tree(covering_vertex_set, rho, size_of_list):
    y = dict()
    n = len(list(covering_vertex_set))
    i = 0
    for j in range(1, n + 1):
        y[j] = covering_vertex_set[i]
        i = i + 1
    pseudo_vertex_set = []
    pseudo_vertex_set.append(covering_vertex_set[0])
    temp_edges_of_spanning_list = []
    counter_2 = 0
    while counter_2 == 0:
        for i in pseudo_vertex_set:
            for j in rho:
                if j != ",":
                    size_of_pseudo_vertex_set = len(pseudo_vertex_set)
                    counter = 0
                    for k in pseudo_vertex_set:
                        if np.array_equal(i.dot(rho.get(j)),k):
                            counter = 1
                    if counter == 0:
                        pseudo_vertex_set.append(i.dot(rho.get(j)))
                        x = i.dot(rho.get(j))
                        for k in range(1, size_of_list + 1):
                            if np.array_equal(y.get(k), i):
                                starting_vertex = k
                            if np.array_equal(y.get(k), x):
                                ending_vertex = k
                        temp_edges_of_spanning_list.append((starting_vertex, ending_vertex))
        len_of_pseudo_vertex_set = len(pseudo_vertex_set)
        len_of_vertex_set = len(covering_vertex_set)
        if len_of_vertex_set == len_of_pseudo_vertex_set:
            counter_2 = 1
    return temp_edges_of_spanning_list


#This function gets us all the paths from base node 1. As we are dealing with a spanning tree, there will be an unique path per end point
def get_all_simple_paths_from_1(spanning_tree_edges, size_of_list):
   temp_vertex_set = []
   for i in range(1, size_of_list+1):
       temp_vertex_set.append(i)
   #print(temp_vertex_set)
   set_of_graph_edges = spanning_tree_edges
   G = nx.DiGraph()
   G.add_nodes_from(temp_vertex_set)
   G.add_edges_from(set_of_graph_edges)
   #print(G)
   all_paths = []
   for i in range(1,size_of_list+1):
       for paths in nx.all_simple_paths(G, source=1, target=i):
           all_paths.append(paths)
   return list(all_paths)

def paths_in_spanning_tree_to_words(list_of_all_paths_from_1, dict_from_loops_to_words, rho, generators_for_free_group):
   list_of_words = []

   #Generate words from loops
   for i in list_of_all_paths_from_1:
       new_word = []
       x = len(i)
       for j in range(x-1):
           a = i[j]
           b = i[j+1]
           edge = (a,b)
           if edge in dict_from_loops_to_words.keys():
               y = dict_from_loops_to_words[edge]
               for key in rho:
                   if np.array_equal(rho[key],y):
                       new_word.append(key)
           else:
               edge = (b,a)
               y = dict_from_loops_to_words[edge]
               for key in rho:
                   if np.array_equal(rho[key], y):
                       list_of_generators = list(generators_for_free_group)
                       for k in list_of_generators:
                           if k != ',':
                               if np.array_equal(key,k):
                                   c = "-"+k
                                   new_word.append(c) # c denote the edge a going in the negative direction
       list_of_words.append(new_word)

   return list_of_words

def construct_loops_in_covering_graph_and_their_words(spanning_tree_path_as_words, list_of_all_paths_from_1, size_of_list, covering_space_edges_set, dict_from_loops_to_words, rho, generators_for_free_group):
    edges_in_spanning_tree = []
    for i in list_of_all_paths_from_1:
        if len(i) == 1:
            edges_in_spanning_tree.append(i)
        else:
            j = 0
            x = len(i)
            for k in range(x-1):
                edges_in_spanning_tree.append([i[k],i[k+1]])

    #Now we remove repeated edges in the set of spanning tree
    edges_in_spanning_tree_without_repeated_element = []
    for i in edges_in_spanning_tree:
        if i not in edges_in_spanning_tree_without_repeated_element:
            edges_in_spanning_tree_without_repeated_element.append(i)

    list_of_path_and_their_inverses = list_of_all_paths_from_1
    list_of_words_and_their_inverses = spanning_tree_path_as_words
    inverse_path = []
    for i in list_of_all_paths_from_1:
        k = i
        temp_inverse_path = k[::-1]
        inverse_path.append(temp_inverse_path)
    inverse_words = []
    for i in spanning_tree_path_as_words:
        #print(i)
        if len(i) == 1:
            temp_word_1 = []
            for k in i:
                temp_word_1.append("-" + str(k))
            inverse_words.append(temp_word_1)
        else:
            temp_word = []
            for k in i[::-1]:
                temp_word.append("-" + str(k))
            inverse_words.append(temp_word)

    inverse_dict = dict() #create dictionary from inverse paths and their words
    for i in inverse_path:
        x = inverse_path.index(i)
        inverse_dict[tuple(i)] = inverse_words[x]

    spanning_tree_path_to_words_dict = dict()
    for i in list_of_all_paths_from_1:
        x = list_of_all_paths_from_1.index(i)
        spanning_tree_path_to_words_dict[tuple(i)] = spanning_tree_path_as_words[x]

    loops_in_covering_graph_word_form = []

    if len(covering_space_edges_set) == len(set(covering_space_edges_set)): #This means that all the edges are distinct
        edges_in_covering_graph_not_spanning_tree_distinct_edges_case = covering_space_edges_set
        for i in edges_in_spanning_tree_without_repeated_element:
            edges_in_covering_graph_not_spanning_tree_distinct_edges_case.remove(tuple(i))

        for i in edges_in_covering_graph_not_spanning_tree_distinct_edges_case:
            i_in_list_form = list(i)
            if i_in_list_form[0] == i_in_list_form[1]:
                if i_in_list_form[0] == 1:
                    matrix_for_loop_at_1 = dict_from_loops_to_words[(1,1)]
                    for j in generators_for_free_group:
                        if j != ",":
                            if numpy.array_equal(rho.get(j), matrix_for_loop_at_1):
                                loops_in_covering_graph_word_form.append(j)
                else:
                    loop = (i_in_list_form[0],i_in_list_form[0])
                    matrix_for_loop = dict_from_loops_to_words[loop]
                    for j in generators_for_free_group:
                        if j != ",":
                            if numpy.array_equal(rho.get(j), matrix_for_loop):
                                generator_for_loop = j
                    for j in list_of_all_paths_from_1:
                        if j[-1] == i_in_list_form[0]:
                            path_to_loop = j
                    path_to_loop_inverse = path_to_loop[::-1]
                    word_for_loop = []
                    word_for_loop = list(spanning_tree_path_to_words_dict[tuple(path_to_loop)]+list(generator_for_loop)+list(inverse_dict[tuple(path_to_loop_inverse)]))
                    loops_in_covering_graph_word_form.append(word_for_loop)
            else:
                starting_vertex = i_in_list_form[0]
                ending_vertex = i_in_list_form[1]
                if ending_vertex != 1:
                    matrix_for_edge = dict_from_loops_to_words[i]
                    for j in generators_for_free_group:
                        if j != ",":
                            if numpy.array_equal(rho.get(j), matrix_for_edge):
                                generator_for_edge = j
                    for j in list_of_all_paths_from_1:
                        if j[-1] == starting_vertex:
                            path_to_starting_vertex_of_non_spanning_tree_edge = j
                        if j[-1] == ending_vertex:
                            path_to_ending_vertex_of_non_spanning_tree_edge = j
                    reverse_path_to_ending_vertex = path_to_ending_vertex_of_non_spanning_tree_edge[::-1]
                    word_for_loop = []
                    word_for_loop = list(spanning_tree_path_to_words_dict[tuple(path_to_starting_vertex_of_non_spanning_tree_edge)])+list(generator_for_edge)+list(inverse_dict[tuple(reverse_path_to_ending_vertex)])
                    loops_in_covering_graph_word_form.append(word_for_loop)
                else:
                    for j in list_of_all_paths_from_1:
                        if j[-1] == starting_vertex:
                            path_to_starting_vertex_of_non_spanning_tree_edge = j
                    matrix_for_edge = dict_from_loops_to_words[i]
                    for j in generators_for_free_group:
                        if j != ",":
                            if numpy.array_equal(rho.get(j), matrix_for_edge):
                                generator_for_edge = j
                    word_for_loop = []
                    word_for_loop = list(spanning_tree_path_to_words_dict[tuple(path_to_starting_vertex_of_non_spanning_tree_edge)])+list(generator_for_edge)
                    loops_in_covering_graph_word_form.append(word_for_loop)
        return loops_in_covering_graph_word_form
    else: #This means there are repeated elements in the set
        #Create edge set without repeated elements
        covering_space_edges_set_no_repeated_elements = []
        repeated_element_in_covering_graph_edges = []
        for i in covering_space_edges_set:
            if i not in covering_space_edges_set_no_repeated_elements:
                covering_space_edges_set_no_repeated_elements.append(i)
            else:
                repeated_element_in_covering_graph_edges.append(i)

        counter = 0
        for i in repeated_element_in_covering_graph_edges:
            if list(i) in edges_in_spanning_tree_without_repeated_element:
                counter = 1

        if counter == 0: #Repeated edges are not inside the spanning tree -- TRICKY CASE
            for i in repeated_element_in_covering_graph_edges:
                i_in_list_form = list(i)
                starting_vertex = i_in_list_form[0]
                ending_vertex = i_in_list_form[1]
                if ending_vertex != 1:
                    if starting_vertex == 1:
                        word_for_loop = []
                        for j in list_of_all_paths_from_1:
                            if j[-1] == ending_vertex:
                                x = j[::-1]
                        word_for_loop = list(spanning_tree_path_to_words_dict[tuple(j)]) + list(inverse_dict[tuple(x)])
                        loops_in_covering_graph_word_form.append(word_for_loop)
                    else:
                        for j in list_of_all_paths_from_1:
                            if j[-1] == starting_vertex:
                                path_to_starting_vertex_of_non_spanning_tree_edge = j
                            if j[-1] == ending_vertex:
                                path_to_ending_vertex_of_non_spanning_tree_edge = j
                            reverse_path_to_ending_vertex = path_to_ending_vertex_of_non_spanning_tree_edge[::-1]
                            word_for_loop = []
                            word_for_loop = list(spanning_tree_path_to_words_dict[tuple(path_to_starting_vertex_of_non_spanning_tree_edge)]) + list(a) + list(inverse_dict[tuple(reverse_path_to_ending_vertex)])
                            loops_in_covering_graph_word_form.append(word_for_loop)
                            word_for_loop = []
                            word_for_loop = list(spanning_tree_path_to_words_dict[tuple(path_to_starting_vertex_of_non_spanning_tree_edge)]) + list(b) + list(inverse_dict[tuple(reverse_path_to_ending_vertex)])
                            loops_in_covering_graph_word_form.append(word_for_loop)
                else:
                    for j in list_of_all_paths_from_1:
                        if j[-1] == starting_vertex:
                            path_to_starting_vertex_of_non_spanning_tree_edge = j
                    word_for_loop = []
                    for k in list(spanning_tree_path_to_words_dict[tuple(path_to_starting_vertex_of_non_spanning_tree_edge)]):
                        word_for_loop.append(k)
                    word_for_loop.append('a')
                    loops_in_covering_graph_word_form.append(word_for_loop)
                    word_for_loop = []
                    for k in list(spanning_tree_path_to_words_dict[tuple(path_to_starting_vertex_of_non_spanning_tree_edge)]):
                        word_for_loop.append(k)
                    word_for_loop.append('b')
                    loops_in_covering_graph_word_form.append(word_for_loop)
            edges_in_covering_graph_with_no_spanning_tree_edges = []
            for i in covering_space_edges_set_no_repeated_elements:
                if list(i) not in edges_in_spanning_tree_without_repeated_element:
                    edges_in_covering_graph_with_no_spanning_tree_edges.append(i)
            edges_in_covering_graph_with_no_spanning_tree_edges_and_no_repeated_loops = []
            for i in edges_in_covering_graph_with_no_spanning_tree_edges:
                if i not in repeated_element_in_covering_graph_edges:
                    edges_in_covering_graph_with_no_spanning_tree_edges_and_no_repeated_loops.append(i)
            for i in edges_in_covering_graph_with_no_spanning_tree_edges_and_no_repeated_loops:
                i_in_list_form = list(i)
                starting_vertex = i_in_list_form[0]
                ending_vertex = i_in_list_form[1]
                if ending_vertex != 1:
                    if starting_vertex == 1:
                        word_for_loop = []
                        for j in list_of_all_paths_from_1:
                            if j[-1] == ending_vertex:
                                x = j[::-1]
                        word_for_loop = list(spanning_tree_path_to_words_dict[tuple(j)]) + list(inverse_dict[tuple(x)])
                        loops_in_covering_graph_word_form.append(word_for_loop)
                    else:
                        matrix_for_edge = dict_from_loops_to_words[i]
                        for j in generators_for_free_group:
                            if j != ",":
                                if numpy.array_equal(rho.get(j), matrix_for_edge):
                                    generator_for_edge = j
                        for j in list_of_all_paths_from_1:
                            if j[-1] == starting_vertex:
                                path_to_starting_vertex_of_non_spanning_tree_edge = j
                            if j[-1] == ending_vertex:
                                path_to_ending_vertex_of_non_spanning_tree_edge = j
                        reverse_path_to_ending_vertex = path_to_ending_vertex_of_non_spanning_tree_edge[::-1]
                        word_for_loop = []
                        word_for_loop = list(spanning_tree_path_to_words_dict[tuple(path_to_starting_vertex_of_non_spanning_tree_edge)]) + list(generator_for_edge) + list(inverse_dict[tuple(reverse_path_to_ending_vertex)])
                        loops_in_covering_graph_word_form.append(word_for_loop)
            return loops_in_covering_graph_word_form
        else:
            #print("One of the repeated edges is inside the spanning tree")
            temp_dict = dict_from_loops_to_words
            # Repeated elements generator should become b in temp dictionary
            for i in repeated_element_in_covering_graph_edges:
                temp_dict[i] = rho.get('b')
            edges_in_covering_graph_not_in_spanning_tree = covering_space_edges_set_no_repeated_elements
            for j in edges_in_spanning_tree_without_repeated_element:
                if tuple(j) not in repeated_element_in_covering_graph_edges:
                    edges_in_covering_graph_not_in_spanning_tree.remove(tuple(j))
            for i in edges_in_covering_graph_not_in_spanning_tree:
                i_in_list_form = list(i)
                if i_in_list_form[0] == i_in_list_form[1]:
                    if i_in_list_form[0] == 1:
                        matrix_for_loop_at_1 = dict_from_loops_to_words[(1, 1)]
                        for j in generators_for_free_group:
                            if j != ",":
                                if numpy.array_equal(rho.get(j), matrix_for_loop_at_1):
                                    loops_in_covering_graph_word_form.append(j)
                    else:
                        loop = (i_in_list_form[0], i_in_list_form[0])
                        matrix_for_loop = dict_from_loops_to_words[loop]
                        for j in generators_for_free_group:
                            if j != ",":
                                if numpy.array_equal(rho.get(j), matrix_for_loop):
                                    generator_for_loop = j
                        for j in list_of_all_paths_from_1:
                            if j[-1] == i_in_list_form[0]:
                                path_to_loop = j
                        path_to_loop_inverse = path_to_loop[::-1]
                        word_for_loop = []
                        word_for_loop = list(spanning_tree_path_to_words_dict[tuple(path_to_loop)] + list(generator_for_loop) + list(inverse_dict[tuple(path_to_loop_inverse)]))
                        loops_in_covering_graph_word_form.append(word_for_loop)
                else:
                    starting_vertex = i_in_list_form[0]
                    ending_vertex = i_in_list_form[1]
                    if ending_vertex != 1:
                        if starting_vertex ==1:
                            word_for_loop = []
                            for j in list_of_all_paths_from_1:
                                if j[-1] == ending_vertex:
                                    x = j[::-1]
                            for k in generators_for_free_group:
                                if k != ",":
                                    if numpy.array_equal(rho.get(k), dict_from_loops_to_words[tuple(i)]):
                                        generator_for_loop = k
                            word_for_loop = list(k)+list(inverse_dict[tuple(x)])
                            loops_in_covering_graph_word_form.append(word_for_loop)
                        else:
                            matrix_for_edge = dict_from_loops_to_words[i]
                            for j in generators_for_free_group:
                                if j != ",":
                                    if numpy.array_equal(rho.get(j), matrix_for_edge):
                                        generator_for_edge = j
                            for j in list_of_all_paths_from_1:
                                if j[-1] == starting_vertex:
                                    path_to_starting_vertex_of_non_spanning_tree_edge = j
                                if j[-1] == ending_vertex:
                                    path_to_ending_vertex_of_non_spanning_tree_edge = j
                            reverse_path_to_ending_vertex = path_to_ending_vertex_of_non_spanning_tree_edge[::-1]
                            word_for_loop = []
                            word_for_loop = list(spanning_tree_path_to_words_dict[tuple(path_to_starting_vertex_of_non_spanning_tree_edge)]) + list(generator_for_edge) + list(inverse_dict[tuple(reverse_path_to_ending_vertex)])
                            loops_in_covering_graph_word_form.append(word_for_loop)
                    else:
                        for j in list_of_all_paths_from_1:
                            if j[-1] == starting_vertex:
                                path_to_starting_vertex_of_non_spanning_tree_edge = j
                        matrix_for_edge = dict_from_loops_to_words[i]
                        for j in generators_for_free_group:
                            if j != ",":
                                if numpy.array_equal(rho.get(j), matrix_for_edge):
                                    generator_for_edge = j
                        word_for_loop = []
                        word_for_loop = list(spanning_tree_path_to_words_dict[tuple(path_to_starting_vertex_of_non_spanning_tree_edge)]) + list(generator_for_edge)
                        loops_in_covering_graph_word_form.append(word_for_loop)
            return loops_in_covering_graph_word_form






        #We use temp dict to get the words for generators of covering graph

def determine_connectedness(covering_graph_visualisation):
    return networkx.is_weakly_connected(covering_graph_visualisation)

def sequence_generator(size_of_sequence):
    x = []
    for i in range(0, size_of_sequence):
        x.append(input("Please enter +1 or -1: \n"))
    return x


def where_fundamental_group_is_sent_to(loops_in_covering_graph_as_words, pseudo_anosov_function):
    x = []
    for i in loops_in_covering_graph_as_words:
        y = []
        if i == str(a):
            x.append(pseudo_anosov_function[a])
        else:
            for j in i:
                if j == str(a):
                    y = y +  pseudo_anosov_function[a]
                if j == str(b):
                    y = y +  pseudo_anosov_function[b]
                inverse_a = "-"+str(a)
                inverse_b = "-" + str(b)
                if j == inverse_a:
                    y = y + pseudo_anosov_function[inverse_a]
                if j == inverse_b:
                    y = y + pseudo_anosov_function[inverse_b]
            x.append(y)
            #print(x)
    #print(x)
    return x



def test_for_list(where_words_in_fundamental_group_is_sent_to, rho, covering_vertex_set):
    """First we create a dictionary to say which matrix the generators a,b,-a and -b map to"""
    generator_and_their_corresponding_matrices = dict()
    for i in rho.keys():
        generator_and_their_corresponding_matrices[i] = rho.get(i)
        x = "-" + str(i)
        generator_and_their_corresponding_matrices[x] = numpy.linalg.inv(rho.get(i))
    #print(generator_and_their_corresponding_matrices)

    """Next we check if the function lifts"""
    size_of_array = len(where_words_in_fundamental_group_is_sent_to)
    counter = 0
    inverse_a = "-" + str(a)
    inverse_b = "-" + str(b)
    for i in where_words_in_fundamental_group_is_sent_to:
        if i == str(a):
            base_vector = covering_vertex_set[0]
            if numpy.array_equal(base_vector.dot(generator_and_their_corresponding_matrices[i]), covering_vertex_set[0]):
                counter = counter + 1
        elif i == str(b):
            base_vector = covering_vertex_set[0]
            if numpy.array_equal(base_vector.dot(generator_and_their_corresponding_matrices[i]), covering_vertex_set[0]):
                counter = counter + 1
        elif i == str(inverse_a):
            base_vector = covering_vertex_set[0]
            if numpy.array_equal(base_vector.dot(generator_and_their_corresponding_matrices[i]), covering_vertex_set[0]):
                counter = counter + 1
        elif i == str(inverse_b):
            base_vector = covering_vertex_set[0]
            if numpy.array_equal(base_vector.dot(generator_and_their_corresponding_matrices[i]), covering_vertex_set[0]):
                counter = counter + 1
        else:
            base_vector = covering_vertex_set[0]
            for j in i:
                base_vector = base_vector.dot(generator_and_their_corresponding_matrices[j])
            if numpy.array_equal(base_vector,covering_vertex_set[0]):
                counter = counter + 1
    if counter == size_of_array:
        return "This passes the test"
    else:
        return "This fails the test"

"""

def build_homological_matrix(covering_space_edges):
    #I am using Covering_space_edges instead of Covering_space_edge_set as for the graphs, covering space edge set doesn't include edges in the spanning tree
    # Need to figure out why
    n = len(covering_space_edges)
    matrix = np.zeros((n, n))
    return matrix

def complete_homological_matrix(homological_matrix, covering_space_edges_set, spanning_tree_edges, pseudo_anosov_function, rho, dict_from_loops_to_words, list_of_all_paths_from_1, size_of_list, generators_for_free_group):
    count = 0
    for e in generators_for_free_group:
        if e != ",":
            count = count + 1

    no_of_edges_in_covering_space = count * size_of_list

    if len(covering_space_edges_set) == no_of_edges_in_covering_space:
        edges_in_covering_graph = covering_space_edges_set
        edges_in_covering_graph = sorted(edges_in_covering_graph)
    else:
        edges_in_covering_graph = covering_space_edges_set + spanning_tree_edges
        edges_in_covering_graph = sorted(edges_in_covering_graph)

    #print("Edges in the covering graph are: " +  str(edges_in_covering_graph))
    #print("The Pseudo-Anosov function is: " + str(pseudo_anosov_function))
    #print("rho is: " + str(rho))
    #print("dict_from_loops_to_words is: " + str(dict_from_loops_to_words))
    #print("list_of_all_paths_from_1: " + str(list_of_all_paths_from_1))
    #print("size of list: "+ str(size_of_list))
    #print("generators of the free group is: " + str(generators_for_free_group))

    temp_rho = dict()

    for i in rho.keys():
        temp_rho[i] = rho[i]
        #print(i)
        matrix_associated_to_key = rho[i]
        #print(matrix_associated_to_key)
        inverse_of_matrix = np.linalg.inv(matrix_associated_to_key)
        #print(inverse_of_matrix)
        #print(generators_for_free_group[0])
        #print(generators_for_free_group[2])
        inverse = "-" + str(i)
        #print(inverse)
        temp_rho[inverse] = inverse_of_matrix


    #print("temp_rho is: " + str(temp_rho))

    key_list = list(rho.keys())

    for i in edges_in_covering_graph:
        #print("edge chosen now is: " + str(i))
        matrix_associated_to_edge = dict_from_loops_to_words[i]
        #print("matrix_associated_to_edge is: " + str(matrix_associated_to_edge))
        for key in rho:
            if np.array_equal(matrix_associated_to_edge, rho[key]):
                word_associated_to_edge = str(key)
        #print("The word_associated_to_edge is: " + str(word_associated_to_edge))
        if i in spanning_tree_edges and i[0] == 1:
            for key in pseudo_anosov_function:
                if str(key) == word_associated_to_edge:
                    word_we_get_after_pseudo_anosov_function_acts_on_i = pseudo_anosov_function[key]
            #print("word_we_get_after_pseudo_anosov_function_acts_on_i is: " + str(word_we_get_after_pseudo_anosov_function_acts_on_i))
            temp_array = []
            alpha = [0]*size_of_list
            alpha[0] = 1
            base_point = numpy.array(alpha)
            #print("basepoint: " + str(base_point))
            starting_vertex_in_the_path_in_covering_space = base_point
            for k in word_we_get_after_pseudo_anosov_function_acts_on_i:
                edge_path = []
                #print("Key is: " + str(k))
                for key in temp_rho:
                    if str(key) == k:
                        for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                            where_is_1 = int(m)+1
                        #print("where_is_1: " + str(where_is_1))
                        edge_path.append(where_is_1)
                        starting_vertex_in_the_path_in_covering_space = starting_vertex_in_the_path_in_covering_space.dot(temp_rho[key])
                        #print("starting_vertex_in_the_path_in_covering_space" + str(starting_vertex_in_the_path_in_covering_space))
                        for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                            where_is_1 = int(m) + 1
                        #print("where is 1" + str(where_is_1))
                        edge_path.append(where_is_1)
                        #print(k)
                        #print("edge path is: " + str(edge_path))
                        if k not in key_list:
                            #print("K not in key list")
                            m = edge_path[0]
                            n = edge_path[1]
                            correct_edge = [-n,-m]
                            #print(correct_edge)
                            #print(edge_path)
                            edge_path = correct_edge
                        temp_array.append(edge_path)
            row_index = edges_in_covering_graph.index(i)
            #print("row_index is" + str(row_index))
            #print("Temp row is: " + str(temp_array))
            for l in temp_array:
                #print(l)
                #column_index = edges_in_covering_graph.index(tuple(l))
                #homological_matrix[row_index][column_index] = temp_array.count(l)
                first_coordinate = l[0]
                second_coordinate = l[1]
                if first_coordinate > 0:
                    column_index = edges_in_covering_graph.index(tuple(l))
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] + 1
                else:
                    new_first_coordinate = tuple([abs(first_coordinate),abs(second_coordinate)])
                    column_index = edges_in_covering_graph.index(new_first_coordinate)
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] - 1


        else:
            edge_i_starting_vertex = i[0]
            #print("edge_i_starting_vertex:" + str(edge_i_starting_vertex))
            for g in list_of_all_paths_from_1:
                if g[-1] == edge_i_starting_vertex:
                    path_to_edge_i = g
            #print("path_to_edge_i: " + str(path_to_edge_i))
            path_we_need = path_to_edge_i+[i[1]]
            #print("path_we_need:" + str(path_we_need))

            #Now we convert the path to words

            new_word = []
            x = len(path_we_need)
            for j in range(x - 1):
                a = path_we_need[j]
                b = path_we_need[j + 1]
                edge = (a, b)
                if edge in dict_from_loops_to_words.keys():
                    y = dict_from_loops_to_words[edge]
                    for key in rho:
                        if np.array_equal(rho[key], y):
                            new_word.append(key)
            #print("new_word" + str(new_word))

            #Now we see where the word is sent to based on the pseudo-Anosov map

            where_is_path_sent_to_by_pseudo_anosov_map = []

            for t in new_word:
                for key in pseudo_anosov_function:
                    if str(key) == t:
                        for r in pseudo_anosov_function[key]:
                            where_is_path_sent_to_by_pseudo_anosov_map.append(r)
            #print("where_is_path_sent_to_by_pseudo_anosov_map:" + str(where_is_path_sent_to_by_pseudo_anosov_map))

            temp_array = []
            alpha = [0] * size_of_list
            alpha[0] = 1
            base_point = numpy.array(alpha)
            starting_vertex_in_the_path_in_covering_space = base_point
            for k in where_is_path_sent_to_by_pseudo_anosov_map:
                edge_path = []
                for key in temp_rho:
                    if str(key) == k:
                        for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                            where_is_1 = int(m) + 1
                        edge_path.append(where_is_1)
                        starting_vertex_in_the_path_in_covering_space = starting_vertex_in_the_path_in_covering_space.dot(temp_rho[key])
                        for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                            where_is_1 = int(m) + 1
                        edge_path.append(where_is_1)
                        if k not in key_list:
                            #print("K not in key list")
                            m = edge_path[0]
                            n = edge_path[1]
                            correct_edge = [-n,-m]
                            #print(correct_edge)
                            #print(edge_path)
                            edge_path = correct_edge
                        temp_array.append(edge_path)
            row_index = edges_in_covering_graph.index(i)
            for l in temp_array:
                #print(l)
                #column_index = edges_in_covering_graph.index(tuple(l))
                #homological_matrix[row_index][column_index] = temp_array.count(l)
                first_coordinate = l[0]
                second_coordinate = l[1]
                if first_coordinate > 0:
                    column_index = edges_in_covering_graph.index(tuple(l))
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] + 1
                else:
                    new_first_coordinate = tuple([abs(first_coordinate),abs(second_coordinate)])
                    column_index = edges_in_covering_graph.index(new_first_coordinate)
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] - 1
    return homological_matrix

def construct_final_homological_matrix(completed_homological_matrix, spanning_tree_edges,covering_space_edges_set, size_of_list, generators_for_free_group):
    temp_homological_matrix = completed_homological_matrix
    count = 0
    for e in generators_for_free_group:
        if e != ",":
            count = count + 1

    no_of_edges_in_covering_space = count * size_of_list
    if len(covering_space_edges_set) == no_of_edges_in_covering_space:
        edges_in_covering_graph = covering_space_edges_set
        edges_in_covering_graph = sorted(edges_in_covering_graph)
    else:
        edges_in_covering_graph = covering_space_edges_set + spanning_tree_edges
        edges_in_covering_graph = sorted(edges_in_covering_graph)
    temp_array = []
    for i in spanning_tree_edges:
        temp_array.append(edges_in_covering_graph.index(i))
    temp_homological_matrix = np.delete(temp_homological_matrix, temp_array, axis=0)
    temp_homological_matrix = np.delete(temp_homological_matrix, temp_array, axis=1)
    return temp_homological_matrix

"""

def build_homological_matrix(covering_space_edges,where_words_in_fundamental_group_is_sent_to):
    #I am using Covering_space_edges instead of Covering_space_edge_set as for the graphs, covering space edge set doesn't include edges in the spanning tree
    # Need to figure out why
    n = len(covering_space_edges)
    m = len(where_words_in_fundamental_group_is_sent_to)
    matrix = np.zeros((m, n))
    #print(matrix)
    return matrix

def complete_homological_matrix(homological_matrix, covering_space_edges_set, spanning_tree_edges, pseudo_anosov_function, rho, dict_from_loops_to_words, list_of_all_paths_from_1, size_of_list, generators_for_free_group, where_words_in_fundamental_group_is_sent_to):
    count = 0
    for e in generators_for_free_group:
        if e != ",":
            count = count + 1

    no_of_edges_in_covering_space = count * size_of_list

    if len(covering_space_edges_set) == no_of_edges_in_covering_space:
        edges_in_covering_graph = covering_space_edges_set
        #edges_in_covering_graph = sorted(edges_in_covering_graph)
    else:
        edges_in_covering_graph = covering_space_edges_set + spanning_tree_edges
        #edges_in_covering_graph = sorted(edges_in_covering_graph)

    #print(edges_in_covering_graph)

    #print("Edges in the covering graph are: " +  str(edges_in_covering_graph))
    #print("The Pseudo-Anosov function is: " + str(pseudo_anosov_function))
    #print("rho is: " + str(rho))
    #print("dict_from_loops_to_words is: " + str(dict_from_loops_to_words))
    #print("list_of_all_paths_from_1: " + str(list_of_all_paths_from_1))
    #print("size of list: "+ str(size_of_list))
    #print("generators of the free group is: " + str(generators_for_free_group))

    temp_rho = dict()

    for i in rho.keys():
        temp_rho[i] = rho[i]
        #print(i)
        matrix_associated_to_key = rho[i]
        #print(matrix_associated_to_key)
        inverse_of_matrix = np.linalg.inv(matrix_associated_to_key)
        #print(inverse_of_matrix)
        #print(generators_for_free_group[0])
        #print(generators_for_free_group[2])
        inverse = "-" + str(i)
        #print(inverse)
        temp_rho[inverse] = inverse_of_matrix
    """

    #print("temp_rho is: " + str(temp_rho))

    key_list = list(rho.keys())


    for i in edges_in_covering_graph:
        #print("edge chosen now is: " + str(i))
        matrix_associated_to_edge = dict_from_loops_to_words[i]
        #print("matrix_associated_to_edge is: " + str(matrix_associated_to_edge))
        for key in rho:
            if np.array_equal(matrix_associated_to_edge, rho[key]):
                word_associated_to_edge = str(key)
        #print("The word_associated_to_edge is: " + str(word_associated_to_edge))
        if i in spanning_tree_edges and i[0] == 1:
            for key in pseudo_anosov_function:
                if str(key) == word_associated_to_edge:
                    word_we_get_after_pseudo_anosov_function_acts_on_i = pseudo_anosov_function[key]
            #print("word_we_get_after_pseudo_anosov_function_acts_on_i is: " + str(word_we_get_after_pseudo_anosov_function_acts_on_i))
            temp_array = []
            alpha = [0]*size_of_list
            alpha[0] = 1
            base_point = numpy.array(alpha)
            #print("basepoint: " + str(base_point))
            starting_vertex_in_the_path_in_covering_space = base_point
            for k in word_we_get_after_pseudo_anosov_function_acts_on_i:
                edge_path = []
                #print("Key is: " + str(k))
                for key in temp_rho:
                    if str(key) == k:
                        for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                            where_is_1 = int(m)+1
                        #print("where_is_1: " + str(where_is_1))
                        edge_path.append(where_is_1)
                        starting_vertex_in_the_path_in_covering_space = starting_vertex_in_the_path_in_covering_space.dot(temp_rho[key])
                        #print("starting_vertex_in_the_path_in_covering_space" + str(starting_vertex_in_the_path_in_covering_space))
                        for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                            where_is_1 = int(m) + 1
                        #print("where is 1" + str(where_is_1))
                        edge_path.append(where_is_1)
                        #print(k)
                        #print("edge path is: " + str(edge_path))
                        if k not in key_list:
                            #print("K not in key list")
                            m = edge_path[0]
                            n = edge_path[1]
                            correct_edge = [-n,-m]
                            #print(correct_edge)
                            #print(edge_path)
                            edge_path = correct_edge
                        temp_array.append(edge_path)
            row_index = edges_in_covering_graph.index(i)
            #print("row_index is" + str(row_index))
            #print("Temp row is: " + str(temp_array))
            for l in temp_array:
                #print(l)
                #column_index = edges_in_covering_graph.index(tuple(l))
                #homological_matrix[row_index][column_index] = temp_array.count(l)
                first_coordinate = l[0]
                second_coordinate = l[1]
                if first_coordinate > 0:
                    column_index = edges_in_covering_graph.index(tuple(l))
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] + 1
                else:
                    new_first_coordinate = tuple([abs(first_coordinate),abs(second_coordinate)])
                    column_index = edges_in_covering_graph.index(new_first_coordinate)
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] - 1


        else:
            edge_i_starting_vertex = i[0]
            #print("edge_i_starting_vertex:" + str(edge_i_starting_vertex))
            for g in list_of_all_paths_from_1:
                if g[-1] == edge_i_starting_vertex:
                    path_to_edge_i = g
            #print("path_to_edge_i: " + str(path_to_edge_i))
            path_we_need = path_to_edge_i+[i[1]]
            #print("path_we_need:" + str(path_we_need))

            #Now we convert the path to words

            new_word = []
            x = len(path_we_need)
            for j in range(x - 1):
                a = path_we_need[j]
                b = path_we_need[j + 1]
                edge = (a, b)
                if edge in dict_from_loops_to_words.keys():
                    y = dict_from_loops_to_words[edge]
                    for key in rho:
                        if np.array_equal(rho[key], y):
                            new_word.append(key)
            #print("new_word" + str(new_word))

            #Now we see where the word is sent to based on the pseudo-Anosov map

            where_is_path_sent_to_by_pseudo_anosov_map = []

            for t in new_word:
                for key in pseudo_anosov_function:
                    if str(key) == t:
                        for r in pseudo_anosov_function[key]:
                            where_is_path_sent_to_by_pseudo_anosov_map.append(r)
            print("where_is_path_sent_to_by_pseudo_anosov_map:" + str(where_is_path_sent_to_by_pseudo_anosov_map))

            temp_array = []
            alpha = [0] * size_of_list
            alpha[0] = 1
            base_point = numpy.array(alpha)
            starting_vertex_in_the_path_in_covering_space = base_point
            for k in where_is_path_sent_to_by_pseudo_anosov_map:
                edge_path = []
                for key in temp_rho:
                    if str(key) == k:
                        for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                            where_is_1 = int(m) + 1
                        edge_path.append(where_is_1)
                        starting_vertex_in_the_path_in_covering_space = starting_vertex_in_the_path_in_covering_space.dot(temp_rho[key])
                        for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                            where_is_1 = int(m) + 1
                        edge_path.append(where_is_1)
                        if k not in key_list:
                            #print("K not in key list")
                            m = edge_path[0]
                            n = edge_path[1]
                            correct_edge = [-n,-m]
                            #print(correct_edge)
                            #print(edge_path)
                            edge_path = correct_edge
                        temp_array.append(edge_path)
            row_index = edges_in_covering_graph.index(i)
            for l in temp_array:
                #print(l)
                #column_index = edges_in_covering_graph.index(tuple(l))
                #homological_matrix[row_index][column_index] = temp_array.count(l)
                first_coordinate = l[0]
                second_coordinate = l[1]
                if first_coordinate > 0:
                    column_index = edges_in_covering_graph.index(tuple(l))
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] + 1
                else:
                    new_first_coordinate = tuple([abs(first_coordinate),abs(second_coordinate)])
                    column_index = edges_in_covering_graph.index(new_first_coordinate)
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] - 1
    """
    key_list = list(rho.keys())
    for i in where_words_in_fundamental_group_is_sent_to:
        #print("word selected is" + str(i))
        where_is_path_sent_to_by_pseudo_anosov_map = i
        #print("where_is_path_sent_to_by_pseudo_anosov_map:" + str(where_is_path_sent_to_by_pseudo_anosov_map))
        temp_array = []
        alpha = [0] * size_of_list
        alpha[0] = 1
        base_point = numpy.array(alpha)
        starting_vertex_in_the_path_in_covering_space = base_point
        for k in where_is_path_sent_to_by_pseudo_anosov_map:
            #print("K is: " + str(k))
            edge_path = []
            for key in temp_rho:
                if str(key) == k:
                    for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                        where_is_1 = int(m) + 1
                    edge_path.append(where_is_1)
                    starting_vertex_in_the_path_in_covering_space = starting_vertex_in_the_path_in_covering_space.dot(
                        temp_rho[key])
                    for m in np.where(starting_vertex_in_the_path_in_covering_space == 1):
                        where_is_1 = int(m) + 1
                    edge_path.append(where_is_1)
                    if k not in key_list:
                        # print("K not in key list")
                        m = edge_path[0]
                        n = edge_path[1]
                        correct_edge = [-n, -m]
                        # print(correct_edge)
                        # print(edge_path)
                        edge_path = correct_edge
                    temp_array.append(edge_path)
                    #print("edge path is: " + str(edge_path))
        row_index = where_words_in_fundamental_group_is_sent_to.index(i)
        #print("where words in fg is send to: " + str(where_words_in_fundamental_group_is_sent_to))
        #print("row index:" + str(row_index))
        #print("edges in covering graph: " + str(edges_in_covering_graph))
        #print("array of edges coreesponding to word is:" + str(temp_array))
        #print("i" + str(i))
        #print(len(temp_array))
        #print(len(i))
        counter = 0
        for l in temp_array:
            #print("l is: " + str(l))
            # column_index = edges_in_covering_graph.index(tuple(l))
            # homological_matrix[row_index][column_index] = temp_array.count(l)
            first_coordinate = l[0]
            second_coordinate = l[1]
            if first_coordinate > 0:
                if edges_in_covering_graph.count(tuple(l)) == 1:
                    column_index = edges_in_covering_graph.index(tuple(l))
                    #print("Column index is: " + str(column_index))
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] + 1
                else:
                    #print(temp_array)
                    index_of_edge_in_temp_array = temp_array.index(l)
                    #print(index_of_edge_in_temp_array)
                    #print(i)
                    #print("index is: " + str(counter))
                    generator_of_i = i[counter]
                    #print("LOOK HERE" + str(generator_of_i))
                    if generator_of_i == 'a':
                        column_index = edges_in_covering_graph.index(tuple(l))
                        #print("Column index is: " + str(column_index))
                        homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] + 1
                    else:
                        column_index = edges_in_covering_graph.index(tuple(l)) + 1
                        #print("Column index is: " + str(column_index))
                        homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] + 1
                #print("Homological matriz is: " + str(homological_matrix))
            else:
                new_first_coordinate = tuple([abs(first_coordinate), abs(second_coordinate)])
                if edges_in_covering_graph.count(new_first_coordinate) == 1:
                    column_index = edges_in_covering_graph.index(new_first_coordinate)
                    #print("Column index is: " + str(column_index))
                    homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] - 1
                else:
                    #print(i)
                    #print("index is: " + str(counter))
                    generator_of_i = i[counter]
                    #print("LOOK HERE" + str(generator_of_i))
                    if generator_of_i == 'a' or generator_of_i == '-a':
                        column_index = edges_in_covering_graph.index(new_first_coordinate)
                        #print("Column index is: " + str(column_index))
                        homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] - 1
                    else:
                        column_index = edges_in_covering_graph.index(new_first_coordinate) + 1
                        #print("Column index is: " + str(column_index))
                        homological_matrix[row_index][column_index] = homological_matrix[row_index][column_index] - 1
                #print("Homological matriz is: " + str(homological_matrix))
            counter = counter + 1
        #print("Homological matriz is: " + str(homological_matrix))
    #print(homological_matrix)
    return homological_matrix

def construct_final_homological_matrix(completed_homological_matrix, spanning_tree_edges,covering_space_edges_set, size_of_list, generators_for_free_group):
    temp_homological_matrix = completed_homological_matrix
    count = 0
    for e in generators_for_free_group:
        if e != ",":
            count = count + 1

    no_of_edges_in_covering_space = count * size_of_list
    if len(covering_space_edges_set) == no_of_edges_in_covering_space:
        edges_in_covering_graph = covering_space_edges_set
        #edges_in_covering_graph = sorted(edges_in_covering_graph)
    else:
        edges_in_covering_graph = covering_space_edges_set + spanning_tree_edges
        #edges_in_covering_graph = sorted(edges_in_covering_graph)
    temp_array = []
    for i in spanning_tree_edges:
        temp_array.append(edges_in_covering_graph.index(i))
    #temp_homological_matrix = np.delete(temp_homological_matrix, temp_array, axis=0)
    temp_homological_matrix = np.delete(temp_homological_matrix, temp_array, axis=1)
    return temp_homological_matrix

def calculate_eigenvalues(final_homological_matrix):
    vals, vecs = LA.eig(final_homological_matrix)
    return vals

def construct_absolute_value_of_eigenvalues(eigenvalues):
    temp_list = []
    for i in eigenvalues:
        temp_list.append(abs(i))
    return temp_list

def getting_the_max_asol_eigenvalue(absolute_value_of_eigenvalues):
    return max(absolute_value_of_eigenvalues)



def main():
   """Contruction of Finite Type Covering Graphs using Monodromy Representation"""

   """We first construct the codomain of the monodromy representation S_n"""
   size_of_list = int(input("What do you want the size of n in S_n to be\n"))
   s_n = list(permutation_list(size_of_list))
   permutation_matrices = build_permutation_matrices(s_n, size_of_list)
   pairs_of_permutation_matrices = list(itertools.combinations(permutation_matrices,2))

   """Here we are asking the user for the input on the generators of the free group so that the monodromy representation can be built"""
   generators_for_free_group = input("What do you want the of generators for free group to be\n")
   free_group = free_group_generator(generators_for_free_group)

   """In this section of the code, we define the Pseudo-Anosov Function"""
   inverse_of_a = "-" + str(a)
   inverse_of_b = "-" + str(b)
   #pseudo_anosov_function = {a: [str(a), str(inverse_of_b) ,str(a)], b: [str(b),str(inverse_of_a)], inverse_of_a: [str(inverse_of_a), str(b), str(inverse_of_a)], inverse_of_b: [str(a), str(inverse_of_b)]}
   #pseudo_anosov_function = {a: [str(a), str(b), str(a)], b: [str(b), str(a)],inverse_of_a: [str(inverse_of_a), str(inverse_of_b), str(inverse_of_a)], inverse_of_b: [str(inverse_of_a), str(inverse_of_b)]}
   pseudo_anosov_function = {a: [str(a), str(b), str(a)], b: [str(a), str(b)],inverse_of_a: [str(inverse_of_a), str(inverse_of_b), str(inverse_of_a)],inverse_of_b: [str(inverse_of_b), str(inverse_of_a)]}
   print("The Psuedo-Anosov function is (In Penner Construction Form) " + str(pseudo_anosov_function))

   table_header = ["Row Number", "rho", "Covering Graph Edges", "Spanning Tree Edges", "Number of Edges in Covering Graph but not Spanning Tree", "List of paths in spanning tree from base point 1", "Spanning Tree Paths as Words", "Generators_for_Fundamental_Group_of_Covering_Graph", "Test to see if all edges not in spanning tree is accounted for", "Where does the pseudo-anosov map send the generators of the fundamental group of the covering space to?", "Can the function be lifted?", "Homological Matrix", "Homological Matrix without spanning tree edges", "eigenvalues", "Absolute Value of Eigenvalues", "Maximum Absolute Eigenvalue" ]
   table = []
   table.append(table_header)
   table_row_number = 1

   """The Following array is defined to be the collection of all spectral radii, one from each cover"""
   set_of_max_eigenvalue_from_each_covering_space = []

   for i in pairs_of_permutation_matrices:
       column_for_table = []
       matrix_1 = i[0]
       matrix_2 = i[1]
       rho = representation_function(generators_for_free_group, matrix_1, matrix_2)

       covering_vertex_set = covering_space_vertices(size_of_list)
       covering_space_edges = generating_covering_space(generators_for_free_group, covering_vertex_set, rho)

       covering_graph_visualisation, covering_space_edges_set, dict_from_loops_to_words = drawing_covering_space(covering_vertex_set, covering_space_edges,rho , generators_for_free_group, size_of_list)
       number_of_edges_for_check_later_on = len(covering_space_edges_set)

       if determine_connectedness(covering_graph_visualisation):
           spanning_tree_edges = create_spanning_tree(covering_vertex_set, rho, size_of_list)
           if len(spanning_tree_edges) != len(covering_space_edges)-size_of_list-1:
               print(str(table_row_number) + "th iteration fails for number of edges in spanning tree")
           list_of_all_paths_from_1 = get_all_simple_paths_from_1(spanning_tree_edges, size_of_list)

           spanning_tree_path_as_words = paths_in_spanning_tree_to_words(list_of_all_paths_from_1, dict_from_loops_to_words, rho, generators_for_free_group)

           loops_in_covering_graph_as_words = construct_loops_in_covering_graph_and_their_words(spanning_tree_path_as_words, list_of_all_paths_from_1, size_of_list, covering_space_edges_set, dict_from_loops_to_words, rho, generators_for_free_group)

           where_words_in_fundamental_group_is_sent_to = where_fundamental_group_is_sent_to(loops_in_covering_graph_as_words, pseudo_anosov_function)

           does_pseudo_anosov_lift = test_for_list(where_words_in_fundamental_group_is_sent_to, rho, covering_vertex_set)

           if does_pseudo_anosov_lift == "This passes the test":
               #print(where_words_in_fundamental_group_is_sent_to)
               column_for_table.append(table_row_number)
               column_for_table.append(rho)
               column_for_table.append(covering_space_edges_set)
               column_for_table.append(spanning_tree_edges)
               column_for_table.append(len(spanning_tree_edges))
               column_for_table.append(list_of_all_paths_from_1)
               column_for_table.append(spanning_tree_path_as_words)
               column_for_table.append(loops_in_covering_graph_as_words)
               if len(loops_in_covering_graph_as_words) == number_of_edges_for_check_later_on - len(spanning_tree_edges):
                   column_for_table.append("Yes")
               column_for_table.append(where_words_in_fundamental_group_is_sent_to)
               column_for_table.append(does_pseudo_anosov_lift)
               homological_matrix = build_homological_matrix(covering_space_edges,where_words_in_fundamental_group_is_sent_to)
               #print(table_row_number)
               #print(homological_matrix)
               completed_homological_matrix = complete_homological_matrix(homological_matrix, covering_space_edges_set,spanning_tree_edges, pseudo_anosov_function,rho, dict_from_loops_to_words,list_of_all_paths_from_1, size_of_list, generators_for_free_group, where_words_in_fundamental_group_is_sent_to)
               #print(completed_homological_matrix)
               column_for_table.append(completed_homological_matrix)
               final_homological_matrix = construct_final_homological_matrix(completed_homological_matrix, spanning_tree_edges,covering_space_edges_set, size_of_list, generators_for_free_group)
               column_for_table.append(final_homological_matrix)
               eigenvalues = calculate_eigenvalues(final_homological_matrix)
               column_for_table.append(eigenvalues)
               absolute_value_of_eigenvalues = construct_absolute_value_of_eigenvalues(eigenvalues)
               column_for_table.append(absolute_value_of_eigenvalues)
               max_absol_eigenvalue = getting_the_max_asol_eigenvalue(absolute_value_of_eigenvalues)
               if max_absol_eigenvalue > 2.65:
                   print(column_for_table)
               set_of_max_eigenvalue_from_each_covering_space.append(max_absol_eigenvalue)




           table.append(column_for_table)
           table_row_number = table_row_number + 1
       else:
           table_row_number = table_row_number + 1

   final_table = tabulate(table, headers='firstrow')
   #print(final_table)
   #print("Total number of covering spaces is: " + str(table_row_number))
   #print("The set consisting of spectral eigenvalues, one from each finite cover: " + str(set_of_max_eigenvalue_from_each_covering_space))
   D_h_f = max(set_of_max_eigenvalue_from_each_covering_space)
   #print("D_h(f) is: " + str(D_h_f))
   #print("The logarithm of D_h(f) is: " + str(math.log(D_h_f)))
   x_axis = set_of_max_eigenvalue_from_each_covering_space
   #print(x_axis)
   y_axis = []
   for i in range(1,len(set_of_max_eigenvalue_from_each_covering_space)+1):
       y_axis.append(i)
   #print(y_axis)
   plt.plot(y_axis,x_axis, "ob")
   plt.xlabel("Finite Covers")
   plt.ylabel("Spectral Radius from Each Finite Type Cover")
   #plt.show()




main()
