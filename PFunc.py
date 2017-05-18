#!/usr/bin/env python3
# sudo apt-get install python3-tk

# This file is part of PFunc. PFunc provides a set of simple tools for users
# to analyze preference functions and other function-valued traits.
#
# Copyright 2016, 2017 Joseph Kilmer
#
# PFunc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PFunc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Import statements
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import tkinter.font as tkFont

from sys import argv
from sys import platform
from os import getcwd
from os import environ
from os import listdir
from os import path
from math import log10
from math import ceil as ceiling
import shelve

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt  # must come after matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvas
from matplotlib.figure import Figure
from datetime import datetime

# If using matplotlib 2+, make it look like matplotlib 1.5.x
if int(matplotlib.__version__.split('.')[0]) >= 2:
    matplotlib.style.use('classic')

# For opening the PDF help file:
if platform == 'win32':
    from os import startfile
else:
    import subprocess

# For finding R on the system:
try:
    import rpy2.robjects as robjects  # must come after matplotlib or numpy
    environ['R_HOME']
except:
    custom_path = '0'
    if 'PFuncPath.txt' in listdir():
        with open('PFuncPath.txt') as pathfile:
            lines = pathfile.readlines()
            for l in lines:
                if l[0:11] == 'custom_path':
                    custom_path = str(l[12:-1])
                    break
    if custom_path == '0':
        if platform == 'win32' and 'R' in listdir('C:\\Program Files'):
            r_versions = []
            for d in listdir('C:\\Program Files\\R'):
                if d[0:2] == 'R-':
                    r_versions.append(d)
            custom_path = 'C:\\Program Files\\R\\' + r_versions[-1]
        elif platform == 'darwin':
            custom_path = '/Library/Frameworks/R.framework/Resources'
        elif platform == 'linux':
            custom_path = '/usr/bin'
        environ['R_HOME'] = custom_path
        environ['R_USER'] = path.dirname(path.realpath(argv[0]))
        import rpy2.robjects as robjects

r = robjects.r


class PrefFunc():
    '''This is the base-level data structure for the program. Each PrefFunc
    object corresponds to an individual in the dataset. This is called when
    opening a new file and when creating group-level splines.

    As input, it takes a dataframe that originated in R, and the names of
    a bunch of different variables that act as settings for generating splines.
    '''
    def __init__(self, r_data_frame, id_number, smoothing_value, current_sp,
                 sp_lim, sp_min, sp_max,
                 loc_peak, peak_min, peak_max,
                 tol_type, tol_drop, tol_absolute, tol_mode,
                 tol_floor, strength_mode, spline_type='individual'):
        self.smoothing_value = smoothing_value
        self.current_sp = current_sp
        self.sp_lim = sp_lim
        self.sp_min = sp_min
        self.sp_max = sp_max
        self.loc_peak = loc_peak
        self.peak_min = peak_min
        self.peak_max = peak_max
        self.tol_type = tol_type
        self.tol_drop = tol_drop
        self.tol_absolute = tol_absolute
        self.tol_mode = tol_mode
        self.tol_floor = tol_floor
        self.strength_mode = strength_mode
        self.r_data_frame = r_data_frame
        self.id_number = id_number
        self.type = spline_type
        self.sp_status = 'magenta'  # magenta = default, cyan = adjusted
        self.update()
        self.name = r('names(%s)[2]' % self.r_data_frame.r_repr())[0]
        self.data_x = r('curr.func$data.x')
        self.data_y = r('curr.func$data.y')
        self.page = ((self.id_number - 1) // 9) + 1
        self.slot = ((self.id_number - 1) % 9) + 1
        self.background = 'white'
        if self.type == 'group':
            self.constituents = r('mydf')
            self.background = '#ffff99'
            self.name = r('names(%s)[3]' % self.r_data_frame.r_repr())[0]

    def update(self):
        self.generate_spline()
        self.populate_stats()

    def generate_spline(self):
        if self.tol_type.get() == 'relative':
            instance_drop = self.tol_drop.get()
            instance_floor = self.tol_floor.get()
        elif self.tol_type.get() == 'absolute':
            instance_drop = 1
            instance_floor = self.tol_absolute.get()
        if self.loc_peak.get() == 0:
            instance_peak = '1'
        elif self.loc_peak.get() == 1:
            instance_peak = 'c(%s, %s)' % (self.peak_min.get(),
                                           self.peak_max.get())
        if self.sp_status == 'magenta':
            self.reset_sp()
        if self.type == 'group':
            r("ind.data <- %s[2:3]" % self.r_data_frame.r_repr())
        else:
            r("ind.data <- %s" % self.r_data_frame.r_repr())
        r("""curr.func <- PFunc(ind.data, 2, %s, peak.within = %s,
                                drop = %s, tol.mode = '%s',
                                sp.binding = %d, min.sp = %s, max.sp = %s,
                                graph.se = TRUE,
                                forgui = TRUE, tol.floor = %s
             )""" % (self.smoothing_value.get(),
                     instance_peak, instance_drop, self.tol_mode.get(),
                     self.sp_lim.get(), self.sp_min.get(), self.sp_max.get(),
                     instance_floor))
        r("master.gam.list[[%s]] <- curr.func$gam.object" % self.id_number)

    def populate_stats(self):
        self.spline_x = r('curr.func$stimulus')
        self.spline_y = r('curr.func$response')
        self.se = r('curr.func$se')
        self.peak_pref = ('%s' % r('curr.func$peak.preference')).split()[1]
        self.peak_resp = ('%s' % r('curr.func$peak.response')).split()[1]
        self.broad_tolerance = ('%s' % r('curr.func$broad.tol')).split()[1]
        self.strict_tolerance = ('%s' % r('curr.func$strict.tol')).split()[1]
        self.broad_tolerance_points = r('curr.func$broad.tol.points')
        self.strict_tolerance_points = r('curr.func$strict.tol.points')
        self.tolerance_height = ('%s' % r('curr.func$tol.height')).split()[1]
        self.hd_strength = ('%s' % r('curr.func$hd.strength')).split()[1]
        self.hi_strength = ('%s' % r('curr.func$hi.strength')).split()[1]
        self.responsiveness = ('%s' % r('curr.func$responsiveness')).split()[1]
        self.axes_ranges = r('range.bundle')  # min.x, max.x, min.y, max.y
        self.smoothing_value.set((
            '%s' % r('curr.func$smoothing.parameter')).split()[1])
        self.is_flat = r('curr.func$is.flat')

    def stiffen(self):
        '''Increase the smoothing parameter'''
        self.smoothing_value.set(self.increment_sp(by=0.1))
        self.sp_status = 'cyan'
        self.update()
        self.current_sp.set(self.smoothing_value.get())

    def loosen(self):
        '''Decrease the smoothing parameter'''
        self.smoothing_value.set(self.increment_sp(by=-0.1))
        self.sp_status = 'cyan'
        self.update()
        self.current_sp.set(self.smoothing_value.get())

    def reset_sp(self):
        '''Reset the smoothing parameter to the default value'''
        self.smoothing_value.set('-1')
        self.sp_status = 'none'  # Protection against infinite loops in update
        self.update()
        self.sp_status = 'magenta'

    def increment_sp(self, by):
        '''Adjust the smoothing parameter by one step up or down.
        Steps are logarithmic.
        '''
        current_sp = float(self.current_sp.get())
        log_sp_val = log10(current_sp)
        round_log_sp_val = round(log_sp_val, 1)
        new_sp_val = round(10 ** (round_log_sp_val + by), 6)
        return str(new_sp_val)

    def update_peak(self):
        '''Update just the peak of the preference function, without running the
        whole PFunc function in R again.
        '''
        previous_peak = self.peak_pref
        if self.loc_peak.get() == 0:
            instance_peak = '1'
        elif self.loc_peak.get() == 1:
            instance_peak = 'c(%s, %s)' % (self.peak_min.get(),
                                           self.peak_max.get())
        peak_bundle = r('''Peak(input.stimuli = %s,
                                preference.function = master.gam.list[[%s]],
                                peak.within = %s,
                                is.flat = %s)
                        ''' % (self.data_x.r_repr(),
                               self.id_number,
                               instance_peak,
                               self.is_flat.r_repr()))
        self.peak_pref = ('%s' % r('%s$peak.preference'
                                   % peak_bundle.r_repr())).split()[1]
        self.peak_resp = ('%s' % r('%s$peak.response'
                                   % peak_bundle.r_repr())).split()[1]
        if self.tol_mode.get() == 'strict' and previous_peak != self.peak_pref:
            self.update_tolerance()

    def update_tolerance(self):
        '''Update just the tolerance of the preference function, without
        running the whole PFunc function in R again.
        '''
        if self.tol_type.get() == 'relative':
            instance_drop = self.tol_drop.get()
            instance_floor = self.tol_floor.get()
        elif self.tol_type.get() == 'absolute':
            instance_drop = 1
            instance_floor = self.tol_absolute.get()
        r('''temp.stim.values <- data.frame(stimulus = %s)
             temp.peak.bundle <- list(peak.preference = %s,
                                      peak.response = %s,
                                      predicting.stimuli = temp.stim.values,
                                      predicted.response = as.vector(%s),
                                      max.stim = max(temp.stim.values),
                                      min.stim = min(temp.stim.values))
          ''' % (self.spline_x.r_repr(),
                 self.peak_pref,
                 self.peak_resp,
                 self.spline_y.r_repr()
                )
        )
        tolerance_bundle = r('''Tolerance(drop = %s,
                                          peak.bundle = temp.peak.bundle,
                                          is.flat = %s,
                                          preference.function =
                                            master.gam.list[[%s]],
                                          tol.floor = %s)
                             ''' % (instance_drop,
                                    self.is_flat.r_repr(),
                                    self.id_number,
                                    instance_floor))
        self.broad_tolerance = ('%s' % r('%s$broad.tolerance'
                                   % tolerance_bundle.r_repr())).split()[1]
        self.strict_tolerance = ('%s' % r('%s$strict.tolerance'
                                   % tolerance_bundle.r_repr())).split()[1]
        self.broad_tolerance_points = r('%s$cross.points'
                                   % tolerance_bundle.r_repr())
        self.strict_tolerance_points = r('%s$strict.points'
                                   % tolerance_bundle.r_repr())
        self.tolerance_height = ('%s' % r('%s$tolerance.height'
                                   % tolerance_bundle.r_repr())).split()[1]


class GraphArea(Frame):
    '''Contains everything in the main viewing window of PFunc, including
    the welcome screen and the graphs.

    Input is particular pieces of display data as well as the names of
    variables controlled by View settings.
    '''
    def __init__(self, individual_dict, current_col, current_page,
                 view_names, view_pts, view_pandtol, view_spline, view_se,
                 tol_mode, input_font, parent=None, **kw):
        Frame.__init__(self, parent, relief=SUNKEN, bd=1)
        self.current_col = current_col
        self.recent_col = IntVar()
        self.current_page = current_page
        self.view_names = view_names
        self.view_pts = view_pts
        self.view_pandtol = view_pandtol
        self.view_spline = view_spline
        self.view_se = view_se
        self.tol_mode = tol_mode
        self.input_font = input_font
        self.parent = parent
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        self.individual_dict = individual_dict
        self.page_dict = {}
        self.individual_slot_dict = {}
        self.slot_dict = {}
        self.tcid_tols = {}
        self.wrapper = Frame(self)
        self.wrapper.grid(row=0, column=0, sticky=NSEW)
        self.wrapper.columnconfigure(0, weight=1)
        self.wrapper.rowconfigure(0, weight=1)
        self.create_welcome()
        self.build_page_controls()
        self.fig = ''
        self.fig_canvas = ''
        self.cid = ''
        self.current_slot = ''
        self.recent_slot = ''
        self.num_pages = 0

    def create_welcome(self):
        self.welcome_canvas = Canvas(self.wrapper, height=550, width=550,
                                     bg='white')
        self.Welcome_text = self.welcome_canvas.create_text(
            275, 225, text='Welcome to PFunc', font=('Helvetica', 36))
        self.instruction_text = self.welcome_canvas.create_text(
            275, 275, text='Open a data file to begin.',
            font=('Helvetica', 12))
        self.copyright_text1 = 'Copyright (C) 2016, 2017 Joseph Kilmer'
        self.copyright_text2 = ('PFunc is distributed under the GNU General '
                                'Public License v3. See About in the Help '
                                'menu for a summary of GPLv3.\nTo view the '
                                'full license, see the accompanying file '
                                'called COPYING.txt or visit '
                                'http://www.gnu.org/licenses/.')
        if platform == 'darwin':
            self.copyright1 = self.welcome_canvas.create_text(
                275, 500, text=self.copyright_text1,
                font=('Helvetica', 10), justify=CENTER)
            self.copyright2 = self.welcome_canvas.create_text(
                275, 530, text=self.copyright_text2,
                font=('Helvetica', 9), justify=CENTER)
        else:
            self.copyright1 = self.welcome_canvas.create_text(
                275, 500, text=self.copyright_text1,
                font=('Helvetica', 8), justify=CENTER)
            self.copyright2 = self.welcome_canvas.create_text(
                275, 530, text=self.copyright_text2,
                font=('Helvetica', 7), justify=CENTER)
        self.welcome_canvas.grid(row=0, column=0, sticky=NSEW)
        self.view = 'welcome'

    def loading_screen(self):
        if self.view == 'welcome':
            self.welcome_canvas.destroy()
        else:
            self.wrapper.destroy()
            self.wrapper = Frame(self)
            self.wrapper.grid(row=0, column=0, sticky=NSEW)
            self.wrapper.columnconfigure(0, weight=1)
            self.wrapper.rowconfigure(0, weight=1)
        self.loading_canvas = Canvas(self.wrapper, height=550, width=550,
                                     bg='gray75')
        self.loading_text = self.loading_canvas.create_text(
            275, 225, text='Loading...', font=('Helvetica', 24))
        self.loading_text2 = self.loading_canvas.create_text(
            275, 260, text='This may take several seconds.', font=('Helvetica', 12))
        self.loading_canvas.lift(self.loading_text)
        self.loading_canvas.grid(row=0, column=0, sticky=NSEW)
        self.loading_canvas.update_idletasks()
        self.view = 'loading'

    def mini_graphs(self, page, and_deselect=True):
        '''Display 3x3 grid of preference function graphs for a given page.'''
        try:
            self.parent.config(cursor='wait')
        except:
            self.parent.config(cursor='watch')
        self.parent.update()
        self.view = 'mini'
        self.slot_dict.clear()
        self.tcid_tols.clear()
        self.individual_slot_dict.clear()
        self.wrapper.destroy()
        self.wrapper = Frame(self)
        self.wrapper.grid(row=0, column=0, sticky=NSEW)
        self.wrapper.columnconfigure(0, weight=1)
        self.wrapper.rowconfigure(0, weight=1)
        if self.first_page_butt.cget('state') == DISABLED:
            self.first_page_butt.configure(state=NORMAL)
            self.back_page_butt.configure(state=NORMAL)
            self.page_num_ent.configure(state=NORMAL)
            self.page_total.configure(state=NORMAL)
            self.next_page_butt.configure(state=NORMAL)
            self.last_page_butt.configure(state=NORMAL)
        self.fig = Figure(figsize=(7, 7))
        self.fig.subplots_adjust(top=0.95, right=0.95, bottom=0.12, hspace=0.4,
                                 wspace=0.3)
        self.fig_canvas = FigureCanvas(self.fig, master=self)
        self.cid = self.fig_canvas.mpl_connect('button_press_event',
                                               self.mini_graph_click)
        self.fig_canvas.get_tk_widget().grid(row=0, column=0, sticky=NSEW)
        # This is what creates the individual graphs:
        counter = 1
        for i in self.page_dict[page]:
            individual = self.individual_dict[i]
            if int(matplotlib.__version__.split('.')[0]) >= 2:
                self.slot_dict[counter] = self.fig.add_subplot(
                    '33%d' % counter, facecolor=individual.background)
            else:
                self.slot_dict[counter] = self.fig.add_subplot(
                    '33%d' % counter, axisbg=individual.background)
            slot = self.slot_dict[counter]
            slot.tick_params(labelsize=10, top='off', right='off')
            slot.spines['top'].set_visible(False)
            slot.spines['right'].set_visible(False)
            self.tcid_tols[str(slot.axes)] = counter
            self.individual_slot_dict[counter] = i
            self.draw_graph(slot, individual)
            counter += 1
        self.fig.text(0.05, 0.45, 'Preference', ha='center', va='bottom',
                      rotation='vertical', fontsize=20)
        self.fig.text(0.53, 0.02, 'Stimulus', ha='center', va='bottom',
                      fontsize=20)
        if self.current_slot != '':
            self.select_mini_graph(self.current_slot, and_deselect)
        self.fig.canvas.draw()
        self.parent.config(cursor='')

    def mega_graph(self, column):
        '''Draw one big graph for a particular individual.'''
        try:
            self.parent.config(cursor='wait')
        except:
            self.parent.config(cursor='watch')
        self.parent.update()
        self.view = 'mega'
        self.wrapper.destroy()
        self.wrapper = Frame(self)
        self.wrapper.grid(row=0, column=0, sticky=NSEW)
        self.fig = Figure(figsize=(7, 7), dpi=80)
        self.fig.subplots_adjust(top=0.95, right=0.95, bottom=0.15, left=0.15,
                                 hspace=0.3, wspace=0.3)
        self.fig_canvas = FigureCanvas(self.fig, master=self)
        self.cid = self.fig_canvas.mpl_connect('button_press_event',
                                               self.mega_graph_click)
        self.fig_canvas.get_tk_widget().grid(row=0, column=0, sticky=NSEW)
        individual = self.individual_dict[column]
        if int(matplotlib.__version__.split('.')[0]) >= 2:
            slot = self.fig.add_subplot('111', facecolor=individual.background)
        else:
            slot = self.fig.add_subplot('111', axisbg=individual.background)
        slot.tick_params(labelsize=20, top='off', right='off', pad=8)
        slot.spines['top'].set_visible(False)
        slot.spines['right'].set_visible(False)
        self.current_page.set(individual.page)
        self.draw_graph(slot, individual)
        self.fig.text(0.05, 0.45, 'Preference', ha='center', va='bottom',
                      rotation='vertical', fontsize=20)
        self.fig.text(0.53, 0.02, 'Stimulus', ha='center', va='bottom',
                      fontsize=20)
        self.fig.canvas.draw()
        self.parent.config(cursor='')

    def mini_graph_click(self, event):
        '''Defines what happens when a mini graph is clicked.
        A single click either selects or deselects the graph.
        A double-click expands the mini graph into a mega graph.
        '''
        if str(event.inaxes) == 'None':
            self.deselect_mini_graph()
            self.current_slot = ''
        elif event.button == 1:
            new_slot = self.tcid_tols[str(event.inaxes)]
            if event.dblclick:
                if self.current_slot == '':
                    self.current_col.set(self.recent_col.get())
                    self.current_slot = self.recent_slot
                self.select_mini_graph(self.current_slot, and_deselect=False)
                self.mega_graph(self.current_col.get())
                self.page_total.configure(text='/ %s' %
                                          len(self.individual_dict))
                self.page_num_ent.configure(textvariable=self.current_col)
            elif self.current_slot != new_slot:
                self.select_mini_graph(new_slot)
                self.current_slot = new_slot
                self.recent_col.set(self.current_col.get())
                self.recent_slot = new_slot
            else:
                self.recent_slot = self.current_slot
                self.recent_col.set(self.current_col.get())
                self.deselect_mini_graph()
                self.current_slot = ''
        self.fig.canvas.draw()

    def mega_graph_click(self, event):
        '''When a mega graph is double-clicked, the view returns to the 3x3
        grid of mini graphs.
        '''
        if event.button == 1 and event.dblclick:
            self.mini_graphs(self.current_page.get(), and_deselect=False)
            self.fig.canvas.draw()
            self.page_total.configure(text='/ %s' % self.num_pages)
            self.page_num_ent.configure(textvariable=self.current_page)

    def select_mini_graph(self, new_slot, and_deselect=True):
        '''Draws a box around a mini graph and displays its stats when the
        mini graph is selected.
        '''
        if and_deselect:
            self.deselect_mini_graph()
        if new_slot != '':
            self.slot_dict[new_slot].spines['bottom'].set_linewidth(2.0)
            self.slot_dict[new_slot].spines['left'].set_linewidth(2.0)
            self.slot_dict[new_slot].spines['top'].set_linewidth(2.0)
            self.slot_dict[new_slot].spines['right'].set_linewidth(2.0)
            self.slot_dict[new_slot].spines['top'].set_visible(True)
            self.slot_dict[new_slot].spines['right'].set_visible(True)
            self.current_col.set(self.individual_slot_dict[new_slot])
        self.event_generate('<<update_sp>>')
        self.event_generate('<<update_summary>>')

    def deselect_mini_graph(self):
        '''Removes the box around the graph and clears the stat display when
        a mini graph is deselected.
        '''
        if self.current_slot != '':
            self.slot_dict[
                self.current_slot].spines['bottom'].set_linewidth(1.0)
            self.slot_dict[self.current_slot].spines['left'].set_linewidth(1.0)
            self.slot_dict[self.current_slot].spines['top'].set_visible(False)
            self.slot_dict[
                self.current_slot].spines['right'].set_visible(False)
            self.current_slot = ''
            self.current_col.set(0)
            self.event_generate('<<clear_display>>')

    def update_graph(self):
        try:
            self.parent.config(cursor='wait')
        except:
            self.parent.config(cursor='watch')
        self.parent.update()
        if self.view == 'mini':
            self.update_mini_graph()
        elif self.view == 'mega':
            self.update_mega_graph()
        self.parent.config(cursor='')

    def update_mini_graph(self):
        '''Draws a new graph in response to changes in settings or smoothing
        parameter.
        '''
        slot = self.current_slot
        if slot != '':
            slot_item = self.slot_dict[slot]
            slot_item.clear()
            slot_item.tick_params(labelsize=10, top='off', right='off')
            individual = self.individual_dict[self.current_col.get()]
            self.tcid_tols[str(slot_item.axes)] = slot
            self.individual_slot_dict[slot] = self.current_col.get()
            self.draw_graph(slot_item, individual)
            self.fig.canvas.draw()

    def update_mega_graph(self):
        '''Draws a new graph in response to changes in settings or smoothing
        parameter.
        '''
        self.fig.clf()
        individual = self.individual_dict[self.current_col.get()]
        if int(matplotlib.__version__.split('.')[0]) >= 2:
            slot = self.fig.add_subplot('111', facecolor=individual.background)
        else:
            slot = self.fig.add_subplot('111', axisbg=individual.background)
        slot.tick_params(labelsize=20, top='off', right='off', pad=8)
        slot.spines['top'].set_visible(False)
        slot.spines['right'].set_visible(False)
        self.draw_graph(slot, individual)
        self.fig.text(0.05, 0.45, 'Preference', ha='center', va='bottom',
                      rotation='vertical', fontsize=20)
        self.fig.text(0.53, 0.02, 'Stimulus', ha='center', va='bottom',
                      fontsize=20)
        self.fig.canvas.draw()

    def draw_graph(self, slot, individual):
        '''Draw a single graph, either in the mini view or the mega view.'''
        slot.axis(list(individual.axes_ranges))
        if self.view == 'mini':
            pt_size = 5
            plt.setp(slot.xaxis.get_majorticklabels(), rotation=60)
        elif self.view == 'mega':
            pt_size = 10
        if self.view_pts.get() == 1 and individual.type == 'individual':
            slot.plot(individual.data_x, individual.data_y, 'k.',
                      markersize=pt_size)
        elif self.view_pts.get() == 1 and individual.type == 'group':
            n_constit = int(r("""
                              tempdf <- %s
                              length(levels(as.factor(tempdf$names)))
                              """ % individual.r_data_frame.r_repr())[0])
            for i in range(0, n_constit):
                r("current.subset.name <- levels(as.factor(tempdf$names))[%d]"
                  % (i + 1))
                r("""current.subset.rows <- which(tempdf$names
                                                  == current.subset.name)""")
                constx = r("tempdf[current.subset.rows, 2]").r_repr()[2: -1]
                constx = eval('[' + constx + ']')
                consty = r("tempdf[current.subset.rows, 3]").r_repr()[2: -1]
                consty = eval('[' + consty + ']')
                slot.plot(constx, consty, color='#cc99ff', linestyle='solid')
        if self.view_pandtol.get() == 1:
            if individual.peak_pref != 'NA':
                slot.plot([individual.peak_pref, individual.peak_pref],
                          [individual.axes_ranges[2], individual.peak_resp],
                          'r-')
            if self.tol_mode.get() == 'broad':
                current_tolerance_points = individual.broad_tolerance_points
            elif self.tol_mode.get() == 'strict':
                current_tolerance_points = individual.strict_tolerance_points
            for t in range(0, len(current_tolerance_points), 2):
                slot.plot([current_tolerance_points[t],
                           current_tolerance_points[t+1]],
                          [individual.tolerance_height,
                           individual.tolerance_height], 'b-')
            # for t in range(0, len(individual.tolerance_points), 2):
            #     slot.plot([individual.tolerance_points[t],
            #                individual.tolerance_points[t+1]],
            #               [individual.tolerance_height,
            #                individual.tolerance_height], 'b-')
        if self.view_spline.get() == 1:
            slot.plot(individual.spline_x, individual.spline_y, 'k-')
        if self.view_se.get() == 1:
            upper_se = []
            lower_se = []
            for i in range(len(individual.se)):
                upper_se.append(individual.spline_y[i] + individual.se[i])
                lower_se.append(individual.spline_y[i] - individual.se[i])
            slot.plot(individual.spline_x, upper_se, color='#666666',
                      linestyle='dashed')
            slot.plot(individual.spline_x, lower_se, color='#666666',
                      linestyle='dashed')
        if self.view_names.get() == 1 and self.view == 'mini':
            slot.set_title(individual.name, size='small')
        elif self.view_names.get() == 1 and self.view == 'mega':
            slot.set_title(individual.name, size='large')
        minx = individual.axes_ranges[0]
        maxx = individual.axes_ranges[1]
        miny = individual.axes_ranges[2]
        maxy = individual.axes_ranges[3]
        dotx = minx - ((maxx - minx) / 10)
        doty = miny - ((maxy - miny) / 10)
        dottype = individual.sp_status
        slot.plot(dotx, doty, color=dottype, marker='.',
                  markersize=(pt_size*2), clip_on=False)
        slot.plot(dotx, doty, color='black', marker='o', fillstyle='none',
                  markersize=(pt_size), clip_on=False)

    def build_page_controls(self):
        '''Initialize the nav buttons at the bottom of the display area.'''
        self.page_controls = Frame(self)
        self.page_controls.grid(row=1, column=0, sticky=EW+S)
        self.page_controls.columnconfigure(0, weight=1)
        self.page_controls.columnconfigure(7, weight=1)
        if platform == 'darwin':
            pd = [8, 8, 0]  # padx first&last, padx back&next, pady for all
        else:
            pd = [1, 4, 0]

        self.first_page_butt = Button(self.page_controls,
                                      text='|<<', padx=pd[0],
                                      pady=pd[2], state=DISABLED,
                                      command=self.first_page)
        self.first_page_butt.grid(row=0, column=1)
        self.back_page_butt = Button(self.page_controls, text='<', padx=pd[1],
                                     pady=pd[2], state=DISABLED,
                                     command=self.back_page)
        self.back_page_butt.grid(row=0, column=2)
        self.page_num_ent = Entry(self.page_controls, width=3, justify=RIGHT,
                                  textvariable=self.current_page,
                                  font=self.input_font)
        self.page_num_ent.grid(row=0, column=3)
        self.page_num_ent.bind('<Return>', self.enter_page_number)
        self.page_num_ent.configure(state=DISABLED)
        self.page_total = Label(self.page_controls, text='/ 0', state=DISABLED)
        self.page_total.grid(row=0, column=4)
        self.next_page_butt = Button(self.page_controls, text='>', padx=pd[1],
                                     pady=pd[2], state=DISABLED,
                                     command=self.next_page)
        self.next_page_butt.grid(row=0, column=5)
        self.last_page_butt = Button(self.page_controls, text='>>|',
                                     padx=pd[0], pady=pd[2], state=DISABLED,
                                     command=self.last_page)
        self.last_page_butt.grid(row=0, column=6)

    def first_page(self):
        '''Jump to the first page'''
        if self.view == 'mini' and self.current_page.get() > 1:
            self.current_page.set(1)
            self.current_slot = ''
            self.mini_graphs(self.current_page.get())
            self.event_generate('<<clear_display>>')
        elif self.view == 'mega' and self.current_col.get() > 1:
            self.current_col.set(1)
            self.recent_col.set(1)
            self.mega_graph(self.current_col.get())
            individual = self.individual_dict[self.current_col.get()]
            self.current_page.set(individual.page)
            self.current_slot = individual.slot
            self.recent_slot = individual.slot
            self.event_generate('<<update_summary>>')
            self.event_generate('<<update_sp>>')

    def back_page(self):
        '''Go back one page'''
        if self.view == 'mini' and self.current_page.get() > 1:
            self.current_page.set(self.current_page.get() - 1)
            self.current_slot = ''
            self.mini_graphs(self.current_page.get())
            self.event_generate('<<clear_display>>')
        elif self.view == 'mega' and self.current_col.get() > 1:
            self.current_col.set(self.current_col.get() - 1)
            self.recent_col.set(self.current_col.get() - 1)
            self.mega_graph(self.current_col.get())
            individual = self.individual_dict[self.current_col.get()]
            self.current_page.set(individual.page)
            self.current_slot = individual.slot
            self.recent_slot = individual.slot
            self.event_generate('<<update_summary>>')
            self.event_generate('<<update_sp>>')

    def next_page(self):
        '''Go forward one page'''
        num_pages = len(self.page_dict)
        num_ind = len(self.individual_dict)
        if self.view == 'mini' and self.current_page.get() < num_pages:
            self.current_page.set(self.current_page.get() + 1)
            self.current_slot = ''
            self.mini_graphs(self.current_page.get())
            self.event_generate('<<clear_display>>')
        elif self.view == 'mega' and self.current_col.get() < num_ind:
            self.current_col.set(self.current_col.get() + 1)
            self.recent_col.set(self.current_col.get() + 1)
            self.mega_graph(self.current_col.get())
            individual = self.individual_dict[self.current_col.get()]
            self.current_page.set(individual.page)
            self.current_slot = individual.slot
            self.recent_slot = individual.slot
            self.event_generate('<<update_summary>>')
            self.event_generate('<<update_sp>>')

    def last_page(self):
        '''Jump ahead to the last page'''
        num_pages = len(self.page_dict)
        num_ind = len(self.individual_dict)
        if self.view == 'mini' and self.current_page.get() < num_pages:
            self.current_page.set(num_pages)
            self.current_slot = ''
            self.mini_graphs(self.current_page.get())
            self.event_generate('<<clear_display>>')
        elif self.view == 'mega' and self.current_col.get() < num_ind:
            self.current_col.set(num_ind)
            self.recent_col.set(num_ind)
            self.mega_graph(self.current_col.get())
            individual = self.individual_dict[self.current_col.get()]
            self.current_page.set(individual.page)
            self.current_slot = individual.slot
            self.recent_slot = individual.slot
            self.event_generate('<<update_summary>>')
            self.event_generate('<<update_sp>>')

    def enter_page_number(self, event):
        '''Executed when the widget is active and the Return key is pressed.
        The graph view updates to the new page specified in the text box.
        '''
        if self.view == 'mini':
            if self.current_page.get() < 1:
                self.current_page.set(1)
            if self.current_page.get() > self.num_pages:
                self.current_page.set(self.num_pages)
            self.deselect_mini_graph()
            self.mini_graphs(self.current_page.get())
        elif self.view == 'mega':
            if self.current_col.get() < 1:
                self.current_col.set(1)
            if self.current_col.get() > len(self.individual_dict):
                self.current_col.set(len(self.individual_dict))
            self.mega_graph(self.current_col.get())
            self.event_generate('<<update_summary>>')
            self.event_generate('<<update_sp>>')


class SmoothingBox(LabelFrame):
    '''The frame for displaying and controlling the smoothing parameter.'''
    def __init__(self, parent=None, text='Smoothing', padx=2, pady=2,
                 heading_font='TkDefaultFont', input_font='TkDefaultFont',
                 row=0, column=0, current_sp='', platform=platform,
                 **kw):
        LabelFrame.__init__(self, parent, text=text,
                            padx=padx, pady=pady,
                            font=heading_font)
        self.grid(row=row, column=column, sticky=EW)
        self.columnconfigure(3, weight=1)

        self.sp_ent = Entry(self, width=10, textvariable=current_sp,
                            state=DISABLED, font=input_font)
        self.sp_ent.grid(row=0, column=0, sticky=W)
        if platform == 'win32':
            self.sp_dn = Button(self, text='-', width=2,
                                command=self.loosen_event,
                                pady=0, state=DISABLED)
            self.sp_up = Button(self, text='+',
                                command=self.stiffen_event,
                                width=2, pady=0, state=DISABLED)
            self.reset_butt = Button(self, text='reset', pady=0,
                                     command=self.reset_sp_event,
                                     state=DISABLED, width=8)
        elif platform == 'darwin':
            self.sp_dn = Button(self, text='-', command=self.loosen_event,
                                padx=8, pady=0, state=DISABLED)
            self.sp_up = Button(self, text='+', command=self.stiffen_event,
                                padx=8, pady=0, state=DISABLED)
            self.reset_butt = Button(self, text='reset', pady=0,
                                     command=self.reset_sp_event,
                                     state=DISABLED)
        else:
            self.sp_dn = Button(self, text='-', command=self.loosen_event,
                                padx=4, pady=0, state=DISABLED)
            self.sp_up = Button(self, text='+', command=self.stiffen_event,
                                padx=2, pady=0, state=DISABLED)
            self.reset_butt = Button(self, text='reset', pady=0,
                                     command=self.reset_sp_event,
                                     state=DISABLED, padx=4)
        self.sp_dn.grid(row=0, column=1, sticky=E)
        self.sp_up.grid(row=0, column=2)
        self.reset_butt.grid(row=0, column=3, sticky=NSEW)
        self.sp_ent.bind('<Return>', self.enter_sp)

    def loosen_event(self):
        self.event_generate('<<loosen>>')

    def stiffen_event(self):
        self.event_generate('<<stiffen>>')

    def reset_sp_event(self):
        self.event_generate('<<reset_sp>>')

    def activate(self):
        self.sp_ent.configure(state=NORMAL)
        self.sp_up.configure(state=NORMAL)
        self.sp_dn.configure(state=NORMAL)
        self.reset_butt.configure(state=NORMAL)

    def enter_sp(self, event):
        self.event_generate('<<enter_sp>>')


class SummaryBox(LabelFrame):
    '''The frame that displays the spline summary for the currently
    selected individual.
    '''
    def __init__(self, parent=None, text='Summary', padx=2, pady=2,
                 heading_font='TkDefaultFont', summary_font='TkDefaultFont',
                 row=0, column=0, **kw):
        LabelFrame.__init__(self, parent, text=text,
                            padx=padx, pady=pady,
                            font=heading_font)
        self.grid(row=row, column=column, sticky=EW)
        self.columnconfigure(0, weight=1)
        self.summary_text = ('Peak Preference:  \n'
                             'Peak Height:  \n'
                             'Tolerance:  \n'
                             'Strength:  \n'
                             'Responsiveness:  \n'
                             'Smoothing:  ')
        self.summary_window = Text(self, height=6, width=25,
                                   font=summary_font)
        self.summary_window.grid(row=0, column=0, sticky=EW)
        self.summary_window.insert(END, self.summary_text)
        self.summary_window.configure(state=DISABLED)

    def update_summary(self, individual=None,
                       strength_mode=None, tol_mode=None):
        self.summary_window.configure(state=NORMAL)
        self.summary_window.delete(1.17, '1.end')
        self.summary_window.delete(2.13, '2.end')
        self.summary_window.delete(3.11, '3.end')
        self.summary_window.delete('4.10', '4.end')
        self.summary_window.delete(5.16, '5.end')
        self.summary_window.delete(6.11, '6.end')
        if individual is not None:
            self.summary_window.insert(1.17, individual.peak_pref)
            self.summary_window.insert(2.13, individual.peak_resp)
            # self.summary_window.insert(3.11, individual.tolerance)
            if tol_mode.get() == 'broad':
                self.summary_window.insert('3.11', individual.broad_tolerance)
            elif tol_mode.get() == 'strict':
                self.summary_window.insert('3.11', individual.strict_tolerance)

            if strength_mode.get() == 'Height-Dependent':
                self.summary_window.insert('4.10', individual.hd_strength)
            elif strength_mode.get() == 'Height-Independent':
                self.summary_window.insert('4.10', individual.hi_strength)
            self.summary_window.insert(5.16, individual.responsiveness)
            self.summary_window.insert(6.11, individual.smoothing_value.get())
        self.summary_window.configure(state=DISABLED)


class ViewBoxItem(Frame):
    '''A single entry in the View Settings frame.'''
    def __init__(self, parent=None, text='', variable=None,
                 pady=0, row=0, column=0, **kw):
        Frame.__init__(self, parent, padx=0, pady=0)
        self.grid(row=row, column=column, sticky=EW)
        self.v_box = Checkbutton(
            self, variable=variable, pady=pady,
            command=lambda: self.event_generate('<<update_all_graphs>>'))
        self.v_box.grid(row=0, column=0)
        self.v_label = Label(self, text=text, pady=pady)
        self.v_label.grid(row=0, column=1, sticky=W)


class ViewBox(LabelFrame):
    '''The frame that contains settings for toggling particular graphical
    elements in the graphs.
    '''
    def __init__(self, parent=None, text='View', padx=2, pady=0,
                 heading_font='TkDefaultFont',
                 view_names_var=None, view_pts_var=None,
                 view_pandtol_var=None, view_spline_var=None,
                 view_se_var=None, row=0, column=0, **kw):
        LabelFrame.__init__(self, parent, text=text,
                            padx=padx, pady=pady,
                            font=heading_font)
        self.grid(row=row, column=column, sticky=EW)
        self.v_names = ViewBoxItem(self, text='Names', pady=pady,
                                   row=0, variable=view_names_var)
        self.v_datapts = ViewBoxItem(self, text='Data Points', pady=pady,
                                     row=1, variable=view_pts_var)
        self.v_pktol = ViewBoxItem(self, text='Peak & Tolerance', pady=pady,
                                   row=2, variable=view_pandtol_var)
        self.v_splines = ViewBoxItem(self, text='Splines', pady=pady,
                                     row=3, variable=view_spline_var)
        self.v_se = ViewBoxItem(self, text='Standard Error', pady=pady,
                                row=4, variable=view_se_var)


class SmoothingLimitsBox(LabelFrame):
    '''The frame that allows users to control the minimum and maximum values for
    smoothing parameters.
    '''
    def __init__(self, parent=None, text='Smoothing Limits', padx=2, pady=2,
                 heading_font='TkDefaultFont', input_font='TkDefaultFont',
                 row=0, column=0, sp_lim_state=None, sp_min=None,
                 sp_max=None, **kw):
        self.sp_lim_state = sp_lim_state
        self.sp_lim_title_frame = Frame()
        LabelFrame.__init__(self, parent, labelwidget=self.sp_lim_title_frame,
                            padx=padx, pady=pady, font=heading_font)
        self.grid(row=row, column=column, sticky=EW)
        self.sp_lim_name = Label(self.sp_lim_title_frame,
                                 text=text, font=heading_font)
        self.sp_lim_name.grid(row=0, column=0)
        self.sp_lim_box = Checkbutton(self.sp_lim_title_frame,
                                      variable=sp_lim_state,
                                      command=self.sp_lim_toggle)
        self.sp_lim_box.grid(row=0, column=1)
        self.sp_lim_min_lab = Label(self, text='   Min')
        self.sp_lim_min_lab.grid(row=0, column=0, sticky=W)
        self.sp_lim_min_ent = Entry(self, width=4,
                                    textvariable=sp_min, font=input_font)
        self.sp_lim_min_ent.grid(row=0, column=1, sticky=W)
        self.sp_lim_max_lab = Label(self, text='   Max')
        self.sp_lim_max_lab.grid(row=0, column=2, sticky=W)
        self.sp_lim_max_ent = Entry(self, width=4,
                                    textvariable=sp_max, font=input_font)
        self.sp_lim_max_ent.grid(row=0, column=3, sticky=W)
        self.sp_lim_min_ent.bind('<Return>', self.enter_sp_lim)
        self.sp_lim_max_ent.bind('<Return>', self.enter_sp_lim)

    def sp_lim_toggle(self, andupdate=TRUE):
        if self.sp_lim_state.get() == 1:
            self.sp_lim_min_ent.configure(state=NORMAL)
            self.sp_lim_min_lab.configure(state=NORMAL)
            self.sp_lim_max_ent.configure(state=NORMAL)
            self.sp_lim_max_lab.configure(state=NORMAL)
        elif self.sp_lim_state.get() == 0:
            self.sp_lim_min_ent.configure(state=DISABLED)
            self.sp_lim_min_lab.configure(state=DISABLED)
            self.sp_lim_max_ent.configure(state=DISABLED)
            self.sp_lim_max_lab.configure(state=DISABLED)
        if andupdate:
            self.event_generate('<<update_magenta_graphs>>')

    def enter_sp_lim(self, event):
        self.event_generate('<<update_magenta_graphs>>')


class LocalPeakBox(LabelFrame):
    '''The frame that allows users to specify the stimulus range to search
    within for local peaks in splines.
    '''
    def __init__(self, loc_peak_state, peak_min, peak_max,
                 parent=None, text='Find Local Peak', padx=2, pady=2,
                 heading_font='TkDefaultFont', input_font='TkDefaultFont',
                 row=0, column=0, **kw):
        self.loc_peak_state = loc_peak_state
        self.peak_title_frame = Frame()
        LabelFrame.__init__(self, parent, labelwidget=self.peak_title_frame,
                            padx=padx, pady=pady)
        self.grid(row=row, column=column, sticky=EW)
        self.peak_name = Label(self.peak_title_frame, text=text,
                               font=heading_font)
        self.peak_name.grid(row=0, column=0)
        self.peak_box = Checkbutton(self.peak_title_frame,
                                    variable=loc_peak_state,
                                    command=self.loc_peak_toggle)
        self.peak_box.grid(row=0, column=1)
        self.peak_btwn_lab1 = Label(self, text='   Between')
        self.peak_btwn_lab1.grid(row=0, column=0)
        self.peak_btwn_ent1 = Entry(self, width=4,
                                    textvariable=peak_min,
                                    font=input_font, state=DISABLED)
        self.peak_btwn_ent1.grid(row=0, column=1, sticky=W)
        self.peak_btwn_lab2 = Label(self, text='and')
        self.peak_btwn_lab2.grid(row=0, column=2)
        self.peak_btwn_ent2 = Entry(self, width=4,
                                    textvariable=peak_max,
                                    font=input_font, state=DISABLED)
        self.peak_btwn_ent2.grid(row=0, column=3, sticky=W)
        self.peak_btwn_ent1.bind('<Return>', self.enter_peak_btwn)
        self.peak_btwn_ent2.bind('<Return>', self.enter_peak_btwn)

    def loc_peak_toggle(self, andupdate=TRUE):
        if self.loc_peak_state.get() == 1:
            self.peak_btwn_ent1.configure(state=NORMAL)
            self.peak_btwn_ent2.configure(state=NORMAL)
        elif self.loc_peak_state.get() == 0:
            self.peak_btwn_ent1.configure(state=DISABLED)
            self.peak_btwn_ent2.configure(state=DISABLED)
        if andupdate:
            #self.event_generate('<<update_all_graphs>>')
            self.event_generate('<<update_all_peaks>>')
            self.event_generate('<<update_summary>>')

    def enter_peak_btwn(self, event):
        #self.event_generate('<<update_all_graphs>>')
        self.event_generate('<<update_all_peaks>>')
        self.event_generate('<<update_summary>>')


class ToleranceBox(LabelFrame):
    '''The frame containing controls for Tolerance.'''
    def __init__(self, tol_type, tol_drop, tol_floor, tol_absolute, tol_mode,
                 parent=None, text='Tolerance', padx=2, pady=0,
                 heading_font='TkDefaultFont', input_font='TkDefaultFont',
                 row=0, column=0, **kw):
        LabelFrame.__init__(self, parent, text=text, padx=padx, pady=pady,
                            font=heading_font)
        self.grid(row=row, column=column, sticky=EW)
        self.tol_type = tol_type
        self.tol_drop = tol_drop
        self.tol_floor = tol_floor
        self.tol_absolute = tol_absolute
        self.tol_mode = tol_mode
        self.tol_rel_sel = Radiobutton(self,
                                       variable=self.tol_type,
                                       value='relative',
                                       command=self.change_tol_type)
        self.tol_rel_sel.grid(row=0, column=0, sticky=E)
        self.tol_rel_lab = Label(self, text='Drop from peak')
        self.tol_rel_lab.grid(row=0, column=1, sticky=W)
        self.tol_rel_ent = Entry(self, width=5,
                                 textvariable=self.tol_drop, font=input_font)
        self.tol_rel_ent.grid(row=0, column=2, sticky=W)
        self.tol_floor_lab = Label(self, text='Floor')
        self.tol_floor_lab.grid(row=1, column=1, sticky=E)
        self.tol_floor_ent = Entry(self, width=5,
                                   textvariable=self.tol_floor,
                                   font=input_font)
        self.tol_floor_ent.grid(row=1, column=2, sticky=W)
        self.tol_abs_sel = Radiobutton(self,
                                       variable=self.tol_type,
                                       value='absolute',
                                       command=self.change_tol_type)
        self.tol_abs_sel.grid(row=2, column=0, sticky=E)
        self.tol_abs_zone = Frame(self)
        self.tol_abs_zone.grid(row=2, column=1, columnspan=2, sticky=W)
        self.tol_abs_lab = Label(self.tol_abs_zone, text='At set value')
        self.tol_abs_lab.grid(row=0, column=0, sticky=W)
        self.tol_abs_ent = Entry(self.tol_abs_zone, width=5,
                                 textvariable=self.tol_absolute,
                                 font=input_font)
        self.tol_abs_ent.grid(row=0, column=1, sticky=W)
        self.tol_mode_zone = Frame(self)
        self.tol_mode_zone.grid(row=3, column=0, sticky=W, columnspan=3)
        self.tol_mode_lab = Label(self.tol_mode_zone, text='Mode')
        self.tol_mode_lab.grid(row=0, column=0, sticky=E)
        self.tol_mode_broad = Radiobutton(
            self.tol_mode_zone, text='Broad', variable=self.tol_mode,
            value='broad',
            command=self.change_tol_mode)
            # command=lambda: self.event_generate('<<update_all_graphs>>'))
        self.tol_mode_broad.grid(row=0, column=1, sticky=W)
        self.tol_mode_stct = Radiobutton(
            self.tol_mode_zone, text='Strict', variable=tol_mode,
            value='strict',
            command=self.change_tol_mode)
        self.tol_mode_stct.grid(row=0, column=2, sticky=W)
        self.tol_rel_ent.bind('<Return>', self.enter_tol_setting)
        self.tol_floor_ent.bind('<Return>', self.enter_tol_setting)
        self.tol_abs_ent.bind('<Return>', self.enter_tol_setting)
        self.change_tol_type(andupdate=False)

    def change_tol_type(self, andupdate=True):
        if self.tol_type.get() == 'relative':
            self.tol_rel_lab.configure(state=NORMAL)
            self.tol_rel_ent.configure(state=NORMAL)
            self.tol_floor_lab.configure(state=NORMAL)
            self.tol_floor_ent.configure(state=NORMAL)
            self.tol_abs_lab.configure(state=DISABLED)
            self.tol_abs_ent.configure(state=DISABLED)
        elif self.tol_type.get() == 'absolute':
            self.tol_rel_lab.configure(state=DISABLED)
            self.tol_rel_ent.configure(state=DISABLED)
            self.tol_floor_lab.configure(state=DISABLED)
            self.tol_floor_ent.configure(state=DISABLED)
            self.tol_abs_lab.configure(state=NORMAL)
            self.tol_abs_ent.configure(state=NORMAL)
        if andupdate:
            #self.event_generate('<<update_all_graphs>>')
            self.event_generate('<<update_all_tolerances>>')
            self.event_generate('<<update_summary>>')

    def enter_tol_setting(self, event):
        #self.event_generate('<<update_all_graphs>>')
        self.event_generate('<<update_all_tolerances>>')
        self.event_generate('<<update_summary>>')

    def change_tol_mode(self):
        self.event_generate('<<update_all_graphs>>')
        self.event_generate('<<update_summary>>')


class StrengthBox(LabelFrame):
    '''The frame that allows users to change between Strength types.'''
    def __init__(self, strength_mode,
                 parent=None, text='Strength', padx=2, pady=2,
                 heading_font='TkDefaultFont', row=0, column=0, **kw):
        LabelFrame.__init__(self, parent, text=text, padx=padx, pady=pady,
                            font=heading_font)
        self.grid(row=row, column=column, sticky=EW)
        self.columnconfigure(0, weight=1)
        self.strength_options = ('Height-Dependent', 'Height-Independent')
        self.strength_selector = OptionMenu(
            self, strength_mode, *self.strength_options,
            command=self.update_summary_event)
        self.strength_selector.grid(row=0, column=0, sticky=EW)

    def update_summary_event(self, strength_option):
        self.event_generate('<<update_summary>>')


class ControlPanel(Frame):
    '''Control Panel contains all the stat readouts and the adjustable
    settings, including the smoothing parameter box, the summary box, the view
    settings, etc.

    Inputs include variable names for all of the settings.
    '''
    def __init__(self, heading_font, input_font, summary_font, current_sp,
                 view_names, view_pts, view_pandtol, view_spline, view_se,
                 sp_lim, sp_min, sp_max,
                 loc_peak, peak_min, peak_max,
                 tol_type, tol_drop, tol_floor, tol_absolute, tol_mode,
                 strength_mode, parent=None, platform=platform, **kw):
        Frame.__init__(self, parent, relief=SUNKEN, bd=1, padx=7, pady=7)
        self.smoothing_box = SmoothingBox(self, heading_font=heading_font,
                                          input_font=input_font,
                                          current_sp=current_sp, row=0)
        self.summary_box = SummaryBox(self, row=1, heading_font=heading_font,
                                      summary_font=summary_font)
        spacer_text = '--Settings--'
        if platform == 'linux':
            spacer_text = '\n' + spacer_text
        self.set_lab = Label(self, text=spacer_text, pady=0, font=heading_font)
        self.set_lab.grid(row=2, column=0)
        self.view_box = ViewBox(self, row=3,
                                view_names_var=view_names,
                                view_pts_var=view_pts,
                                view_pandtol_var=view_pandtol,
                                view_spline_var=view_spline,
                                view_se_var=view_se,
                                heading_font=heading_font)
        self.smoothing_limits_box = SmoothingLimitsBox(
            self, row=4, heading_font=heading_font, input_font=input_font,
            sp_lim_state=sp_lim, sp_min=sp_min, sp_max=sp_max)
        self.peak_box = LocalPeakBox(parent=self, loc_peak_state=loc_peak,
                                     peak_min=peak_min, peak_max=peak_max,
                                     heading_font=heading_font,
                                     input_font=input_font, row=5)
        self.tolerance_box = ToleranceBox(parent=self, tol_type=tol_type,
                                          tol_drop=tol_drop,
                                          tol_floor=tol_floor,
                                          tol_absolute=tol_absolute,
                                          tol_mode=tol_mode,
                                          heading_font=heading_font,
                                          input_font=input_font,
                                          row=6)
        self.strength_box = StrengthBox(parent=self,
                                        strength_mode=strength_mode,
                                        heading_font=heading_font, row=7)

    def update_summary(self, individual=None,
                       strength_mode=None, tol_mode=None):
        self.summary_box.update_summary(individual, strength_mode, tol_mode)

    def activate(self):
        self.smoothing_box.activate()
        self.active_mode = 'activated'


class FileMenu(Menubutton):
    '''Defines the File menu at the top of the screen (and accompanying
    functions).
    '''
    def __init__(self, file_opt, parent=None, row=0, column=0):
        Menubutton.__init__(self, parent, text='File')
        self.grid(row=row, column=column, sticky=W)
        self.file_opt = file_opt
        self.parent = parent
        self.primary_menu = Menu(self, tearoff=0)
        self.open_menu = Menu(self, tearoff=0)
        self.open_menu.add_command(label='Horizontal...',
                                   command=self.open_horizontal_file)
        self.open_menu.add_command(label='Vertical...',
                                   command=self.open_vertical_file)
        self.primary_menu.add_cascade(label='Open Data File',
                                      menu=self.open_menu)
        self.primary_menu.add_separator()
        self.primary_menu.add_command(label='Load Smoothing Values...',
                                      command=self.open_sp,
                                      state=DISABLED)
        self.primary_menu.add_command(label='Save Smoothing Values...',
                                      command=self.save_sp,
                                      state=DISABLED)
        self.primary_menu.add_command(label='Clear Smoothing Values',
                                      command=self.clear_sps,
                                      state=DISABLED)
        self.primary_menu.add_separator()
        self.primary_menu.add_command(label='Load Previous Settings',
                                      command=self.open_sett)
        self.primary_menu.add_command(label='Save Current Settings',
                                      command=self.save_sett)
        self.primary_menu.add_command(label='Restore Default Settings',
                                      command=self.reset_sett)
        self.primary_menu.add_separator()
        self.primary_menu.add_command(label='Output Spline Figures...',
                                      command=self.output_graphs,
                                      state=DISABLED)
        self.primary_menu.add_command(label='Output Spline Summaries...',
                                      command=self.output_summaries,
                                      state=DISABLED)
        self.primary_menu.add_command(label='Output Spline Points...',
                                      command=self.output_points,
                                      state=DISABLED)
        self.primary_menu.add_command(label='Output Tolerance Points...',
                                      command=self.output_tol,
                                      state=DISABLED)
        self.primary_menu.add_separator()
        self.primary_menu.add_command(label='Quit', command=self.quit)

        self['menu'] = self.primary_menu

    def activate_menu_options(self):
        self.primary_menu.entryconfigure(2, state=NORMAL)
        self.primary_menu.entryconfigure(3, state=NORMAL)
        self.primary_menu.entryconfigure(4, state=NORMAL)
        self.primary_menu.entryconfigure(10, state=NORMAL)
        self.primary_menu.entryconfigure(11, state=NORMAL)
        self.primary_menu.entryconfigure(12, state=NORMAL)
        self.primary_menu.entryconfigure(13, state=NORMAL)

    def _check_missing_stim(self, is_vertical=0):
        '''Used when opening a new file. Checks whether any x-axis values
        are missing.
        '''
        r = robjects.r
        if not is_vertical:
            stim_column = '1'
        else:
            stim_column = 'stim.column'
        missing_stim = int(r('as.numeric(InCheck(NA, mydata[, %s]))'
                             % stim_column)[0])
        if missing_stim:
            error_text = ("Could not open the data file because there "
                          "seems to be one or more missing stimulus "
                          "values.")
            messagebox.showerror('Error', error_text)
            self.event_generate('<<add_message>>', x=106)
            return False
        else:
            return True

    def _check_num_datapoints(self, is_vertical):
        '''Used when opening a new file. Checks whether there are enough data
        points.
        '''
        r = robjects.r
        r("minimum_datapoints <- 10")
        if not is_vertical:
            min_pts = int(r("""for (n in 1:length(name.vect)) {
                                 response <- mydata[, (n + 1)]
                                 num_datapoints = sum(!is.na(response))
                                 minimum_datapoints <- min(minimum_datapoints,
                                                           num_datapoints)
                               }
                               return(minimum_datapoints)
                            """)[0])
        else:
            min_pts = int(r("""for (n in 1:length(name.vect)) {
                                 response <- mydata[, resp.column][which(
                                   mydata[, id.column] == name.vect[n])]
                                 num_datapoints = sum(!is.na(response))
                                 minimum_datapoints <- min(minimum_datapoints,
                                                           num_datapoints)
                               }
                               return(minimum_datapoints)
                            """)[0])
        if min_pts < 3:
            self.event_generate('<<add_message>>', x=103)
            error_text = ("Not enough data to work with. PFunc needs a "
                          "minimum of three data points to make a single "
                          "spline. Make sure that each individual has at "
                          "least three responses.\n\n"
                          "If you are trying to make a preference function "
                          "by combining responses from different individuals, "
                          "then you should group those individuals together "
                          "in your data file. See README file "
                          "(or Help) for more."
            )
            messagebox.showerror('Error', error_text)
            return False
        return True

    def check_data_formatting(self, datafile=None):
        '''Used when opening a new data file. Checks whether file is .csv'''
        if datafile is None:
            return False
        data_formatted_correctly = int(r("""mydata <- read.csv("%s")
                                            if (ncol(mydata) == 1){
                                              mydata <- read.delim("%s")
                                            }
                                            if (ncol(mydata) == 1){
                                              return("0")
                                            } else {
                                              return("1")
                                            }
                                            """ % (datafile.name,
                                                   datafile.name))[0])
        if not data_formatted_correctly:
            error_text = ("The data file you selected is not formatted "
                          "correctly.\n\nMake sure it is saved as a "
                          ".csv file.")
            messagebox.showerror('Error', error_text)
            self.event_generate('<<add_message>>', x=104)
            return False
        else:
            return True

    def open_horizontal_file(self):
        datafile = filedialog.askopenfile(mode='r', **self.file_opt)
        if self.check_data_formatting(datafile):
            r = robjects.r
            r("name.vect = names(mydata)[2: ncol(mydata)]")
            is_vertical = 0
            if (self._check_missing_stim(is_vertical) &
                    self._check_num_datapoints(is_vertical)):
                r("""
                    max.resp <- max(mydata[ , 2:ncol(mydata)], na.rm = TRUE)
                    min.resp <- min(mydata[ , 2:ncol(mydata)], na.rm = TRUE)
                    resp.range <- max.resp - min.resp
                    max.y <- max.resp + (0.0375 * resp.range * 2)
                    min.y <- min.resp - (0.0375 * resp.range * 1)

                    max.stim <- max(mydata[ , 1], na.rm = TRUE)
                    min.stim <- min(mydata[ , 1], na.rm = TRUE)
                    stim.range <- max.stim - min.stim
                    max.x <- max.stim + (0.0375 * stim.range * 1)
                    min.x <- min.stim - (0.0375 * stim.range * 1)

                    range.bundle <- c(min.x, max.x, min.y, max.y)
                    """)
                self.event_generate('<<open_data_file>>', x=is_vertical)

    def open_vertical_file(self):
        datafile = filedialog.askopenfile(mode='r', **self.file_opt)
        if self.check_data_formatting(datafile):
            self.id_column = StringVar()
            self.stim_column = StringVar()
            self.resp_column = StringVar()
            popup = DataDefiner(datafile, self.id_column, self.stim_column,
                                self.resp_column, return_to=self,
                                parent=self.parent.parent)

    def open_vertical_file2(self):
        r = robjects.r
        r("id.column <- which(names(mydata) == '%s')" % self.id_column.get())
        r("stim.column <- which(names(mydata) == '%s')"
            % self.stim_column.get())
        r("resp.column <- which(names(mydata) == '%s')"
            % self.resp_column.get())
        r("""name.vect <- vector()
             for (r in 1:nrow(mydata)){
                ind.id <- as.character(mydata[, id.column][r])
                if (!InCheck(ind.id, name.vect)){
                    name.vect = append(name.vect, ind.id)
                }
             }
          """)
        is_vertical = 1
        if (self._check_missing_stim(is_vertical) &
                self._check_num_datapoints(is_vertical)):
            r("""
                max.resp <- max(mydata[, resp.column])
                min.resp <- min(mydata[, resp.column])
                resp.range <- max.resp - min.resp
                max.y <- max.resp + 0.0375 * resp.range
                min.y <- min.resp - 0.0375 * resp.range

                max.stim <- max(mydata[, stim.column])
                min.stim <- min(mydata[, stim.column])
                stim.range <- max.stim - min.stim
                max.x <- max.stim + 0.0375 * stim.range
                min.x <- min.stim - 0.0375 * stim.range

                range.bundle <- c(min.x, max.x, min.y, max.y)
                """)
            self.event_generate('<<open_data_file>>', x=is_vertical)

    def open_sp(self):
        self.event_generate('<<open_smoothing_file>>')

    def save_sp(self):
        self.event_generate('<<save_smoothing_values>>')

    def clear_sps(self):
        self.event_generate('<<clear_smoothing_values>>')

    def open_sett(self):
        self.event_generate('<<load_settings>>')

    def save_sett(self):
        self.event_generate('<<save_settings>>')

    def reset_sett(self):
        self.event_generate('<<reset_settings>>')

    def output_graphs(self):
        self.event_generate('<<output_graphs>>')

    def output_summaries(self):
        self.event_generate('<<output_summaries>>')

    def output_points(self):
        self.event_generate('<<output_points>>')

    def output_tol(self):
        self.event_generate('<<output_tol>>')

    def quit(self):
        self.event_generate('<<quit>>')


class AdvancedMenu(Menubutton):
    '''Defines the Advanced menu at the top of the screen (and accompanying
    functions).
    '''
    def __init__(self, parent=None, row=0, column=0):
        Menubutton.__init__(self, parent, text='Advanced')
        self.grid(row=row, column=column, sticky=W)
        self.primary_menu = Menu(self, tearoff=0)
        self.primary_menu.add_command(label='Show Message Log',
                                      command=self.message_log)
        self.primary_menu.add_command(label='Construct Group-Level Spline...',
                                      command=self.construct_group_spline,
                                      state=DISABLED)
        self['menu'] = self.primary_menu

    def activate_menu_options(self):
        self.primary_menu.entryconfigure(1, state=NORMAL)

    def message_log(self):
        self.event_generate('<<open_message_log>>')

    def construct_group_spline(self):
        self.event_generate('<<open_group_spline_window>>')


class HelpMenu(Menubutton):
    '''Defines the Help menu at the top of the screen (and accompanying
    functions).
    '''
    def __init__(self, parent=None, row=0, column=0):
        Menubutton.__init__(self, parent, text='Help')
        self.grid(row=row, column=column, sticky=W)
        self.primary_menu = Menu(self, tearoff=0)
        self.primary_menu.add_command(label='Help', command=self.open_help)
        self.primary_menu.add_command(label='About', command=self.about_window)
        self['menu'] = self.primary_menu

    def about_window(self):
        self.event_generate('<<create_about_window>>')

    def open_help(self):
        try:
            if platform == 'win32':
                startfile('README.pdf')
            elif platform == 'darwin':
                subprocess.call(['open', 'README.pdf'])
            else:  # linux
                subprocess.call(['xdg-open', 'README.pdf'])
        except:
            warning_text = ("PFunc failed to locate and open README.pdf. "
                            "You can download a copy of this help file "
                            "from github.com/joccalor/pfunc")
            self.warning = messagebox.showwarning('Warning', warning_text)


class MenuBar(Frame):
    '''Defines the entire menu bar at the top of the screen.'''
    def __init__(self, file_opt, parent=None, row=0, column=0):
        Frame.__init__(self, parent)
        self.parent = parent
        self.grid(row=row, column=column, sticky=EW, columnspan=2)
        self.columnconfigure(3, weight=1)
        self.file_menu = FileMenu(parent=self, file_opt=file_opt)
        self.advc_menu = AdvancedMenu(self, column=1)
        self.help_menu = HelpMenu(self, column=2)

    def activate(self):
        self.file_menu.activate_menu_options()
        self.advc_menu.activate_menu_options()


class PFuncToplevel(Toplevel):
    '''A generic popup window for PFunc (a superclass)'''
    def __init__(self, parent=None, **kw):
        Toplevel.__init__(self, parent, takefocus=True, **kw)
        try:
            img = PhotoImage(file='PFuncIcon.png')
            self.tk.call('wm', 'iconphoto', self._w, img)
        except:
            a = 1


class DataDefiner(PFuncToplevel):
    '''Used in opening a vertical file. Asks users to specify which columns
    of the data contain certain data types.
    '''
    def __init__(self, datafile, id_column, stim_column, resp_column,
                 return_to, parent=None, **kw):
        PFuncToplevel.__init__(self, parent)
        self.return_to = return_to
        self.datafile = datafile
        self.transient(parent)
        rootWd = int(parent.winfo_width()) / 4
        rootHt = int(parent.winfo_height()) / 3
        self.XPos = int(parent.winfo_geometry().split('+')[1]) + rootWd
        self.YPos = int(parent.winfo_geometry().split('+')[2]) + rootHt
        self.geometry('+%d+%d' % (self.XPos, self.YPos))
        self.name_label = Label(self, text='Individual IDs: ')
        self.name_label.grid(row=0, column=0)
        self.xdata_label = Label(self, text='Stimuli (x-axis): ')
        self.xdata_label.grid(row=1, column=0)
        self.ydata_label = Label(self, text='Responses (y-axis): ')
        self.ydata_label.grid(row=2, column=0)
        self.column_names = list(r("names(mydata)"))
        self.id_column = id_column
        self.stim_column = stim_column
        self.resp_column = resp_column
        self.id_column.set('Select')
        self.stim_column.set('Select')
        self.resp_column.set('Select')
        self.column_menu1 = OptionMenu(self, self.id_column,
                                       *self.column_names)
        self.column_menu1.grid(row=0, column=1)
        self.column_menu2 = OptionMenu(self, self.stim_column,
                                       *self.column_names)
        self.column_menu2.grid(row=1, column=1)
        self.column_menu3 = OptionMenu(self, self.resp_column,
                                       *self.column_names)
        self.column_menu3.grid(row=2, column=1)
        self.spacer = Frame(self)
        self.spacer.grid(row=3, column=0, columnspan=2)
        self.okay_butt = Button(self, text='Okay', command=self.okay)
        self.okay_butt.grid(row=4, column=0)
        self.cancel_butt = Button(self, text='Cancel', command=self.cancel)
        self.cancel_butt.grid(row=4, column=1)
        self.column_defs = {'name': self.id_column,
                            'stim': self.stim_column,
                            'resp': self.resp_column}

    def cancel(self):
        self.destroy()

    def okay(self):
        if self.id_column.get() != 'Select' and\
                self.id_column.get() != self.stim_column.get() and\
                self.id_column.get() != self.resp_column.get() and\
                self.stim_column.get() != 'Select' and\
                self.stim_column.get() != self.resp_column.get() and\
                self.resp_column.get() != 'Select':
            self.destroy()
            self.return_to.open_vertical_file2()
        else:
            warning_text = ("You must select columns in your data that "
                            "correspond to each of these three categories.")
            self.warning = messagebox.showwarning('Warning', warning_text)


class GroupSplineWindow(PFuncToplevel):
    '''Used for combining multiple splines into one group-level spline. Users
    tell PFunc which individuals should be part of the group.
    '''
    def __init__(self, parent, individual_dict, combomode, input_font, **kw):
        self.parent = parent
        PFuncToplevel.__init__(self, self.parent)
        self.transient(self.parent)
        self.individual_dict = individual_dict
        self.combomode = combomode
        self.columnconfigure(0, weight=1)
        self.rowconfigure(2, weight=1)
        rootWd = int(parent.winfo_width()) / 2
        rootHt = int(parent.winfo_height()) / 2
        reqWd = int(self.winfo_reqwidth())
        reqHt = int(self.winfo_reqheight())
        XPos = int(parent.winfo_geometry().split('+')[1]) + rootWd - reqWd
        YPos = int(parent.winfo_geometry().split('+')[2]) + rootHt - reqHt
        self.geometry('+%d+%d' % (XPos, YPos))
        self.newname = StringVar()
        self.newname.set('spline%s' % str(len(individual_dict) + 1))
        self.instructions = ("Select the individuals to be used\n"
                             "in this group-level spline.\n\n"
                             "Click and drag to select multiple\n"
                             "individuals at once. Hold down\n"
                             "ctrl to add or subtract individuals\n"
                             "from your selection.")
        self.instruction_box = Label(self, text=self.instructions,
                                     justify=LEFT, padx=5, pady=5)
        self.instruction_box.grid(row=0, column=0, sticky=W)

        self.namebox = Frame(self, pady=10, padx=20)
        self.namebox.grid(row=1, column=0, sticky=EW)
        self.namebox.columnconfigure(1, weight=1)

        self.newname_lab = Label(self.namebox, text='Name')
        self.newname_lab.grid(row=1, column=0, sticky=EW)
        self.newname_ent = Entry(self.namebox, textvariable=self.newname,
                                 width=15, font=input_font)
        self.newname_ent.grid(row=1, column=1, sticky=W)

        self.listframe = Frame(self, padx=20)
        self.listframe.grid(row=2, column=0, sticky=NSEW)
        self.listframe.columnconfigure(0, weight=1)
        self.listframe.rowconfigure(0, weight=1)

        self.namestring = ''
        for i in self.individual_dict:
            self.namestring = (self.namestring +
                               self.individual_dict[i].name + ' ')
        self.namestring = self.namestring[:-1]
        self.names = StringVar()
        self.names.set(self.namestring)
        self.listscroll = Scrollbar(self.listframe, orient=VERTICAL)
        self.listscroll.grid(row=0, column=1, sticky=NS+W)
        self.listbox = Listbox(self.listframe, listvariable=self.names,
                               height=15, selectmode=EXTENDED,
                               yscrollcommand=self.listscroll.set,
                               font=input_font)
        self.listbox.grid(row=0, column=0, sticky=NSEW)
        self.listscroll['command'] = self.listbox.yview

        self.okayframe = Frame(self, padx=20, pady=5)
        self.okayframe.grid(row=3, column=0)
        self.okay_butt = Button(self.okayframe, text='Okay', command=self.okay)
        self.okay_butt.grid(row=0, column=0, sticky=E)
        self.cancel_butt = Button(self.okayframe, text='Cancel',
                                  command=self.cancel)
        self.cancel_butt.grid(row=0, column=1, sticky=W)
        self.event_generate('<<open_message_log>>')

    def cleanup_name(self):
        name = self.newname.get()
        new_name = ""
        alphabet = "abcdefghijklmnopqrstuvwxyz"
        replace_with_underscore = """ <>()[]{}#"'=+-!@#$%^&*`~,\|/?"""
        if name[0].lower() not in alphabet:
            name = 'x' + name[:]
        for c in name:
            if c in replace_with_underscore:
                new_name += '_'
            else:
                new_name += c
        self.newname.set(new_name)

    def cancel(self):
        self.destroy()

    def okay(self):
        self.cleanup_name()
        self.output_dict = {'name': self.newname.get(), 'individual_nums': [],
                            'method': self.combomode.get(), }
        for i in self.listbox.curselection():
            self.output_dict['individual_nums'].append(i+1)
        numsExist = len(self.output_dict['individual_nums']) > 0
        namesExist = len(self.output_dict['name']) > 0
        if numsExist and namesExist:
            self.combine_spline()
            self.destroy()

    def combine_spline(self):
        r('mylist <- list()')
        for i in self.listbox.curselection():
            tempind = self.individual_dict[i + 1]
            tempx = str(tempind.spline_x.r_repr())
            tempy = str(tempind.spline_y.r_repr())
            r("""mylist$%s <- list('xvals' = %s,
                'yvals' = %s)""" % (tempind.name, tempx, tempy))
        if self.combomode.get() == 'none':
            r("""
                xvalues <- vector()
                yvalues <- vector()
                names <- vector()
                for (i in 1:length(mylist)) {
                    xvalues <- c(xvalues, mylist[[i]]$xvals)
                    yvalues <- c(yvalues, mylist[[i]]$yvals)
                    names <- c(names, rep(names(mylist)[i],
                                          length(mylist[[i]]$xvals)))
                }
                mydf <- data.frame(names=names, xvalues=xvalues,
                                   %s=yvalues, stringsAsFactors=FALSE)
            """ % self.newname.get())
        else:
            r("""
                mydf <- data.frame(xvals=NA)
                for(i in mylist){
                  for(j in i$xvals){
                    if(InCheck(j, mydf$xvals)){
                    } else {
                      mydf <- rbind(mydf, j)
                    }
                  }
                }
                for(i in 1:length(mylist)){
                  mydf <- cbind(mydf, NA)
                  names(mydf)[i+1] <- paste('col', i, sep='')
                  for(j in 1:length(mylist[[i]]$xvals)){
                    row <- which(mydf$xvals == mylist[[i]]$xvals[j])
                    mydf[row, i+1] <- mylist[[i]]$yvals[j]
                  }
                }
                mydf <- mydf[2:nrow(mydf), ]

                mydf <- cbind(mydf, NA, NA)
                names(mydf)[length(names(mydf))-1] <- "n"
                names(mydf)[length(names(mydf))] <- "%s"
                for(i in 1:nrow(mydf)){
                  jvec <- vector()
                  for(j in 2:(ncol(mydf)-2)){
                    if(!is.na(mydf[i, j])){
                      jvec[length(jvec) + 1] <- mydf[i, j]
                    }
                    mydf$n[i] <- length(jvec)
                    mydf$%s[i] <- %s(jvec)
                  }
                }
            """ % (self.newname.get(), self.newname.get(),
                   self.combomode.get()))
        self.parent.event_generate('<<add_group_spline>>')


class PFuncMessages(PFuncToplevel):
    '''Defines the popup window of messages that users can access from the
    Advanced menu.
    '''
    def __init__(self, parent, messages, *kw):
        self.parent = parent
        PFuncToplevel.__init__(self, self.parent)
        self.messages = messages
        self.title('PFunc Message Log')
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.logArea = Text(self, height=8, width=32, wrap=WORD)
        self.logArea.insert(END, self.messages.get())
        self.logArea.grid(row=0, column=0, sticky=NSEW)
        self.logArea.tag_add('message_tag', '@0,0', END)
        self.logArea.tag_config('message_tag', lmargin2='32p')
        self.logArea.see(END)
        self.logScroll = Scrollbar(self, orient=VERTICAL,
                                   command=self.logArea.yview)
        self.logScroll.grid(row=0, column=1, sticky=NS+E)
        self.logArea['yscrollcommand'] = self.logScroll.set
        self.logArea.configure(state=DISABLED)
        self._establish_placement()

    def _establish_placement(self):
        screenWd = int(self.parent.winfo_screenwidth())
        reqWd = int(self.winfo_reqwidth())
        reqHt = int(self.winfo_reqheight())
        rootWd = int(self.parent.winfo_width())
        rootHt = int(self.parent.winfo_height()) / 2
        root_leftbound = int(self.parent.winfo_geometry().split('+')[1])
        root_rightbound = int(root_leftbound + rootWd)
        left_space = root_leftbound
        right_space = screenWd - root_rightbound
        if right_space > left_space:
            xOption1 = root_rightbound + 20
            xOption2 = screenWd - reqWd - 20
            xPos = min(xOption1, xOption2)
        else:
            xOption1 = 20
            xOption2 = root_leftbound - reqWd - 50
            xPos = max(xOption1, xOption2)
        yPos = int(self.parent.winfo_geometry().split('+')[2]) + rootHt - reqHt
        self.geometry('+%d+%d' % (xPos, yPos))

    def add_message(self, message_string):
        self.logArea.configure(state=NORMAL)
        self.logArea.insert(END, message_string)
        self.logArea.tag_add('message_tag', '@0,0', END)
        self.logArea.tag_config('message_tag', lmargin2='32p')
        self.logArea.see(END)
        self.logArea.configure(state=DISABLED)


class AboutWindow(PFuncToplevel):
    def __init__(self, parent, title_font, *kw):
        PFuncToplevel.__init__(self, padx=5, pady=5)
        self.parent = parent
        self.title_font = title_font
        self.transient()
        self.title('About PFunc')
        self.initiate_text()
        self.place_elements()
        self.set_geometry()

    def initiate_text(self):
        self.title_text = 'PFunc'
        self.subtitle_text = ('A tool for analyzing preference functions and '
                              'other function-valued traits.\n')
        self.version_text = 'version 0.10.0 \n (2017-05-18)\n'
        self.copyright_text = ('Copyright (C) 2016, 2017 Joseph Kilmer \n\n'
                               'PFunc is distributed under the GNU General '
                               'Public License v3. A full copy of\n'
                               'the license is available in the accompanying '
                               'file called COPYING.txt.\n\n'
                               'PFunc is free software: you can redistribute '
                               'it and/or modify it under the\n'
                               'terms of the GNU General Public License as '
                               'published by the Free Software\n'
                               'Foundation, either version 3 of the License, '
                               'or (at your option) any\n'
                               'later version.\n\n'
                               'PFunc is distributed in the hope that it will '
                               'be useful, but WITHOUT ANY\n'
                               'WARRANTY; without even the implied warranty '
                               'of MERCHANTABILITY or FITNESS FOR\n'
                               'A PARTICULAR PURPOSE. See the GNU General '
                               'Public License for more details.\n\n'
                               'You should have received a copy of the GNU '
                               'General Public License along with\n'
                               'this program. If not, see '
                               'http://www.gnu.org/licenses/.\n')

    def place_elements(self):
        self.title = Label(self, text=self.title_text, font=self.title_font)
        self.title.grid(row=0, column=0)
        self.subtitle = Label(self, text=self.subtitle_text)
        self.subtitle.grid(row=1, column=0)
        self.version = Label(self, text=self.version_text)
        self.version.grid(row=2, column=0)
        self.copyright = Label(self, text=self.copyright_text)
        self.copyright.grid(row=3, column=0)
        self.closebutton = Button(self, text='Close', command=self.destroy)
        self.closebutton.grid(row=4, column=0)

    def set_geometry(self):
        rootWd = int(self.parent.root.winfo_width()) / 2
        rootHt = int(self.parent.root.winfo_height()) / 2
        reqWd = int(self.winfo_reqwidth())
        reqHt = int(self.winfo_reqheight())
        XPos = (int(self.parent.root.winfo_geometry().split('+')[1]) +
                rootWd - reqWd)
        YPos = (int(self.parent.root.winfo_geometry().split('+')[2]) +
                rootHt - reqHt)
        self.geometry('+%d+%d' % (XPos, YPos))


class MainApp():
    '''This is the wrapper for the whole program. It contains and governs
    all the individual pieces.
    '''
    def __init__(self):
        self.root = Tk()
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        self._setup_fonts()
        self._setup_dicts()
        self._setup_variables()
        self._setup_message_lookup()
        self._setup_R()
        self._setup_file_opt()
        self._setup_event_bindings()
        self._setup_window_geometry()
        self.settings_to_default()
        self.menu_bar = MenuBar(file_opt=self.file_opt, parent=self.root)
        self.graph_zone = GraphArea(self.individual_dict, self.current_col,
                                    self.current_page, self.view_names,
                                    self.view_pts, self.view_pandtol,
                                    self.view_spline, self.view_se,
                                    self.tol_mode,
                                    self.input_font, parent=self.root)
        self.graph_zone.grid(row=1, column=0, sticky=NSEW)
        self.control_panel = ControlPanel(heading_font=self.heading_font,
                                          input_font=self.input_font,
                                          summary_font=self.summary_font,
                                          current_sp=self.current_sp,
                                          view_names=self.view_names,
                                          view_pts=self.view_pts,
                                          view_pandtol=self.view_pandtol,
                                          view_spline=self.view_spline,
                                          view_se=self.view_se,
                                          sp_lim=self.sp_lim,
                                          sp_min=self.sp_min,
                                          sp_max=self.sp_max,
                                          loc_peak=self.loc_peak,
                                          peak_min=self.peak_min,
                                          peak_max=self.peak_max,
                                          tol_type=self.tol_type,
                                          tol_drop=self.tol_drop,
                                          tol_floor=self.tol_floor,
                                          tol_absolute=self.tol_absolute,
                                          tol_mode=self.tol_mode,
                                          strength_mode=self.strength_mode,
                                          parent=self.root)
        self.control_panel.grid(row=1, column=1, sticky=NSEW)
        self.root.title('PFunc')
        try:
            img = PhotoImage(file='PFuncIcon.png')
            self.root.tk.call('wm', 'iconphoto', self.root._w, img)
        except:
            a = 1
        self.root.event_generate('<<add_message>>', x=100)
        self.root.config(cursor='')

    def _setup_fonts(self):
        self.default_font = tkFont.nametofont('TkDefaultFont')
        self.small_font = self.default_font.copy()
        self.heading_font = self.default_font.copy()
        self.summary_font = tkFont.nametofont('TkFixedFont')
        self.summary_font.configure(size=9)
        self.input_font = self.default_font.copy()
        self.about_font1 = self.default_font.copy()
        self.small_font.configure(size=8)
        self.heading_font.configure(size=9)
        self.about_font1.configure(size=20)
        if platform == 'win32':
            self.default_font.configure(size=8)
            self.small_font.configure(size=7)
            self.heading_font.configure(size=8)
            self.summary_font = self.default_font.copy()
        elif platform == 'darwin':
            self.default_font.configure(size=10)
            self.small_font.configure(size=9)
            self.heading_font.configure(size=11)
            self.input_font.configure(size=10)
            self.summary_font = self.input_font
            # self.summary_font = tkFont.nametofont('TkTextFont')
            # self.summary_font.configure(size=10)
        else:  # including platform == 'linux'
            self.input_font = tkFont.nametofont('TkTextFont')

    def _setup_dicts(self):
        self.sp_dict = {}  # A dictionary of smoothing parameters
        self.individual_dict = {}  # A dictionary of PrefFunc objects

    def _setup_variables(self):
        self.view_pts = IntVar()
        self.view_pandtol = IntVar()
        self.view_spline = IntVar()
        self.view_names = IntVar()
        self.view_se = IntVar()
        self.sp_lim = IntVar()
        self.sp_min = StringVar()
        self.sp_max = StringVar()
        self.loc_peak = IntVar()
        self.peak_min = StringVar()
        self.peak_max = StringVar()
        self.tol_type = StringVar()
        self.tol_drop = StringVar()
        self.tol_absolute = StringVar()
        self.tol_mode = StringVar()
        self.tol_floor = StringVar()
        self.strength_mode = StringVar()

        self.combomode = StringVar()
        self.messages = StringVar()

        self.current_sp = StringVar()
        self.current_page = IntVar()
        self.current_col = IntVar()
        self.file_type = StringVar()

        self.current_page.set(0)

        self.vertColResp = StringVar()

    def _setup_message_lookup(self):
        self.message_lookup = {}
        self.message_lookup[100] = ("Welcome to PFunc. Open a data file to "
                                    "begin. See the Help menu or the README "
                                    "file for help.")
        self.message_lookup[101] = "Cleared previous smoothing values."
        self.message_lookup[102] = "Opened a new file."
        self.message_lookup[103] = ("Refused to open file because "
                                    "one or more individuals had fewer than "
                                    "three data points.")
        self.message_lookup[104] = ("Failed to open file because it was not "
                                    "properly formatted.")
        self.message_lookup[105] = ("One or more individuals in this "
                                    "dataset have fewer than 10 data points. "
                                    "Consider lowering the minimum smoothing "
                                    "value limit.")
        self.message_lookup[106] = ("Failed to open file because there are "
                                    "fewer stimuli than responses.")

    def _setup_R(self):
        current_directory = StringVar()  # For some reason it must be StrinVar
        current_directory.set(getcwd())
        if platform == 'win32':
            current_directory.set(path.dirname(path.realpath(argv[0])))
            current_directory.set(current_directory.get().replace("\\", "/"))
        r("setwd('%s')" % current_directory.get())
        r("source('PFunc_RCode.R')")

    def _setup_file_opt(self):
        self.file_opt = {}
        self.file_opt['defaultextension'] = '.csv'
        self.file_opt['filetypes'] = [('all files', '.*'),
                                      ('csv files', '.csv'),
                                      ('text files', '.txt')]
        self.file_opt['parent'] = self.root
        self.file_opt['title'] = 'Select a file...'

    def _setup_event_bindings(self):
        self.root.bind('<<open_data_file>>', self.open_data_file)
        self.root.bind('<<update_summary>>', self.update_summary)
        self.root.bind('<<clear_display>>', self.clear_display)
        self.root.bind('<<update_sp>>', self.update_sp)
        self.root.bind('<<loosen>>', self.loosen)
        self.root.bind('<<stiffen>>', self.stiffen)
        self.root.bind('<<reset_sp>>', self.reset_sp)
        self.root.bind('<<enter_sp>>', self.enter_sp)
        self.root.bind('<<update_all_graphs>>', self.update_all_graphs)
        self.root.bind('<<update_all_peaks>>', self.update_all_peaks)
        self.root.bind('<<update_all_tolerances>>', self.update_all_tolerances)
        self.root.bind('<<update_magenta_graphs>>', self.update_magenta_graphs)
        self.root.bind('<<open_message_log>>', self.open_message_log)
        self.root.bind('<<add_message>>', self.add_message)
        self.root.bind('<<open_group_spline_window>>',
                       self.open_group_spline_window)
        self.root.bind('<<add_group_spline>>', self.add_group_spline)
        self.root.bind('<<open_smoothing_file>>', self.open_smoothing_file)
        self.root.bind('<<save_smoothing_values>>', self.save_smoothing_values)
        self.root.bind('<<clear_smoothing_values>>',
                       self.clear_smoothing_values)
        self.root.bind('<<load_settings>>', self.load_settings)
        self.root.bind('<<save_settings>>', self.save_settings)
        self.root.bind('<<reset_settings>>', self.reset_settings)
        self.root.bind('<<output_graphs>>', self.output_graphs)
        self.root.bind('<<output_summaries>>', self.output_summaries)
        self.root.bind('<<output_points>>', self.output_points)
        self.root.bind('<<output_tol>>', self.output_tol)
        self.root.bind('<<quit>>', self.quit)
        self.root.bind('<<create_about_window>>', self.create_about_window)

    def _setup_window_geometry(self):
        if platform == 'darwin':
            scwd = (self.root.winfo_screenwidth() - 717) / 2
            if self.root.winfo_screenheight() > 612:
                scht = (self.root.winfo_screenheight() - 612) / 2
            else:
                scht = self.root.winfo_screenheight()
        elif platform == 'win32':
            scwd = (self.root.winfo_screenwidth() - 746) / 2
            if self.root.winfo_screenheight() > 660:
                scht = (self.root.winfo_screenheight() - 600 - 60) / 2
            else:
                scht = self.root.winfo_screenheight()
        else:
            scwd = (self.root.winfo_screenwidth() - 767) / 2
            if self.root.winfo_screenheight() > 607:
                scht = (self.root.winfo_screenheight() - 607) / 2
            else:
                scht = self.root.winfo_screenheight()
        self.root.geometry('+%d+%d' % (scwd, scht))
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)

    def settings_to_default(self, event=None):
        self.view_names.set(0)
        self.view_pts.set(1)
        self.view_pandtol.set(1)
        self.view_spline.set(1)
        self.view_se.set(0)
        self.tol_type.set('relative')
        self.tol_drop.set('1/3')
        self.tol_absolute.set('1')
        self.tol_mode.set('broad')
        self.tol_floor.set('0')
        self.loc_peak.set(0)
        self.peak_min.set('min')
        self.peak_max.set('max')
        if r("InCheck('min.stim', objects())")[0]:
            self.peak_min.set(r("min.stim")[0])
            self.peak_max.set(r("max.stim")[0])
        self.strength_mode.set('Height-Dependent')
        self.sp_lim.set(1)
        self.sp_min.set('0.05')
        self.sp_max.set('5')
        self.combomode.set('none')

    def _check_num_datapoints(self):
        r = robjects.r
        minimum_datapoints = 10
        for i in self.individual_dict.values():
            r('checkdata <- %s' % i.r_data_frame.r_repr())
            num_datapoints = int(r('sum(!is.na(checkdata[, 2]))')[0])
            minimum_datapoints = min(minimum_datapoints, num_datapoints)
        if minimum_datapoints < 10:
            self.root.event_generate('<<add_message>>', x=105)

    def open_data_file(self, event=None):
        r = robjects.r
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        # self.root.update()
        self.graph_zone.current_slot = ''
        self.graph_zone.page_dict.clear()
        self.graph_zone.individual_dict.clear()
        self.sp_dict.clear()
        r("master.gam.list <- list()")
        if self.graph_zone.view == 'mini' or self.graph_zone.view == 'mega':
            self.root.event_generate('<<add_message>>', x=101)
        self.graph_zone.loading_screen()
        if event.x == 0:
            self.file_type.set('horizontal')
            num_ind = r("ncol(mydata)")[0] - 1
        elif event.x == 1:
            self.file_type.set('vertical')
            num_ind = r("length(name.vect)")[0]
        self.peak_min.set(r("min.stim")[0])
        self.peak_max.set(r("max.stim")[0])
        for i in range(1, num_ind + 1):
            self.sp_dict[i] = StringVar()
            self.sp_dict[i].set('-1')
            if self.file_type.get() == 'horizontal':
                r("""individual_df <- data.frame(stimulus = mydata[, 1],
                    response = mydata[, (%s + 1)])
                    """ % i)
            elif self.file_type.get() == 'vertical':
                r("""individual_df <- data.frame(
                    stimulus = mydata[, stim.column][which(mydata[, id.column]
                                                           == name.vect[%d])],
                    response = mydata[, resp.column][which(mydata[, id.column]
                                                           == name.vect[%d])])
                    """ % (i, i))
            individual_df = r("""
                # individual_df <- data.frame(stimulus = mydata[, 1],
                #                             response = mydata[, (%s + 1)])
                #names(individual_df)[2] <- name.vect[#s]
                individual_name <- name.vect[%s]
                individual_name_char1 <- strsplit(individual_name, "")[[1]][1]
                allowed_characters <-
                  "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ."
                allowed_characters_split <- strsplit(allowed_characters,
                                                     "")[[1]]
                if (!InCheck(individual_name_char1, allowed_characters_split)){
                  individual_name <- paste("X", individual_name, sep="")
                }
                names(individual_df)[2] <- individual_name
                rejector <- vector()
                for (r in 1:nrow(individual_df)) {
                  if (is.na(individual_df[r, 2])){
                    rejector <- c(rejector, r)
                  }
                }
                if (length(rejector) > 0) {
                  individual_df <- individual_df[-rejector, ]
                }
                individual_df
            """ % (i, i))
            self.individual_dict[i] = PrefFunc(
                individual_df, i, self.sp_dict[i], self.current_sp,
                self.sp_lim, self.sp_min, self.sp_max,
                self.loc_peak, self.peak_min, self.peak_max,
                self.tol_type, self.tol_drop, self.tol_absolute, self.tol_mode,
                self.tol_floor, self.strength_mode)
        self.clear_display()
        self.num_pages = num_ind//9
        if num_ind//9 != num_ind/9:
            self.num_pages += 1
        self.graph_zone.num_pages = self.num_pages
        for p in range(1, (self.num_pages + 1)):
            if p < self.num_pages:
                self.graph_zone.page_dict[p] = list(range(1+9*(p-1),
                                                          10+9*(p-1)))
            else:
                remaining_ind = num_ind - (p-1)*9
                ind_list = []
                for r in range(1, (remaining_ind + 1)):
                    ind_list.append(r+9*(p-1))
                self.graph_zone.page_dict[p] = ind_list
        self.graph_zone.mini_graphs(1)
        self.graph_zone.page_total.configure(text='/ %s' % self.num_pages)
        self.current_page.set(1)
        self.control_panel.activate()
        self.menu_bar.activate()
        self.root.event_generate('<<add_message>>', x=102)
        self._check_num_datapoints()
        self.root.config(cursor='')
        # self.root.update()

    def update_summary(self, event=None):
        if self.current_col.get() != 0:
            current_individual = self.individual_dict[self.current_col.get()]
            self.sp_dict[self.current_col.get()] = \
                current_individual.smoothing_value
        else:
            current_individual = None
        self.control_panel.update_summary(individual=current_individual,
                                          strength_mode=self.strength_mode,
                                          tol_mode=self.tol_mode)

    def update_sp(self, event=None):
        if self.current_col.get() != 0:
            self.current_sp.set(self.individual_dict[
                                self.current_col.get()].smoothing_value.get())
        else:
            self.current_sp.set('')

    def clear_display(self, event=None):
        self.current_sp.set('')
        self.current_col.set(0)
        self.update_summary(event=None)

    def loosen(self, event=None):
        col = self.current_col.get()
        if col != 0:
            self.individual_dict[col].loosen()
            self.update_summary(event=None)
            self.graph_zone.update_graph()

    def stiffen(self, event=None):
        col = self.current_col.get()
        if col != 0:
            self.individual_dict[col].stiffen()
            self.update_summary(event=None)
            self.graph_zone.update_graph()

    def reset_sp(self, event=None):
        col = self.current_col.get()
        if col != 0:
            self.individual_dict[col].reset_sp()
            self.graph_zone.update_graph()
            self.update_summary(event=None)
            self.current_sp.set(
                self.individual_dict[col].smoothing_value.get())

    def enter_sp(self, event=None):
        col = self.current_col.get()
        if col != 0:
            self.sp_dict[col].set(self.current_sp.get())
            self.individual_dict[col].sp_status = 'cyan'
            self.individual_dict[col].update()
            self.graph_zone.update_graph()
            self.update_summary(event=None)

    def update_all_graphs(self, event=None):
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        if self.graph_zone.view == 'mini':
            self.graph_zone.mini_graphs(self.current_page.get(),
                                        and_deselect=False)
        elif self.graph_zone.view == 'mega':
            self.graph_zone.mega_graph(self.current_col.get())
        self.root.config(cursor='')

    def update_all_peaks(self, event=None):
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        for i in self.individual_dict:
            self.individual_dict[i].update_peak()
        self.update_all_graphs()
        self.root.config(cursor='')

    def update_all_tolerances(self, event=None):
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        for i in self.individual_dict:
            self.individual_dict[i].update_tolerance()
        self.update_all_graphs()
        self.root.config(cursor='')

    def update_magenta_graphs(self, event=None):
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        if self.graph_zone.num_pages > 0:
            for i in self.individual_dict:
                if self.individual_dict[i].sp_status == 'magenta':
                    # sp_lim_on = (self.sp_lim.get() == 1)
                    # sp_too_small = (
                    #     self.individual_dict[i].smoothing_value.get()
                    #     < self.sp_min.get())
                    # sp_too_big = (
                    #     self.individual_dict[i].smoothing_value.get()
                    #     > self.sp_max.get())
                    # if sp_lim_on and (sp_too_small or sp_too_big):
                    #     self.individual_dict[i].reset_sp()
                    # elif not sp_lim_on:
                    #     self.individual_dict[i].reset_sp()
                    self.individual_dict[i].reset_sp()
            self.update_summary(self.current_col.get())
            if self.graph_zone.view == 'mini' and self.current_col.get() != 0:
                self.graph_zone.select_mini_graph(self.graph_zone.current_slot,
                                                  and_deselect=False)
            self.update_all_graphs()
            self.root.config(cursor='')

    def open_message_log(self, event=None):
        self.logWindow = PFuncMessages(self.root, self.messages)

    def add_message(self, event=None):
        message_code = event.x
        message_string = self.message_lookup[message_code]
        current_datetime = str(datetime.now())
        spc_indx = current_datetime.find(" ")
        time_str = current_datetime[spc_indx + 1: spc_indx+6]
        if self.messages.get() == '':
            message_string = time_str + ' ' + message_string
        else:
            message_string = '\n' + time_str + ' ' + message_string
        self.messages.set(self.messages.get() + message_string)
        for child in self.root.winfo_children():
            if type(child) == PFuncMessages:
                child.add_message(message_string)

    def open_group_spline_window(self, event=None):
        group_spline_window = GroupSplineWindow(self.root,
                                                self.individual_dict,
                                                self.combomode,
                                                self.input_font)

    def add_group_spline(self, event=None):
        newsplinedf = r('mydf')
        self.sp_dict[(len(self.sp_dict) + 1)] = StringVar()
        self.sp_dict[len(self.sp_dict)].set('-1')
        self.individual_dict[(len(self.individual_dict) + 1)] = \
            PrefFunc(newsplinedf, len(self.individual_dict) + 1,
                     self.sp_dict[len(self.sp_dict)], self.current_sp,
                     self.sp_lim, self.sp_min, self.sp_max,
                     self.loc_peak, self.peak_min, self.peak_max,
                     self.tol_type, self.tol_drop, self.tol_absolute,
                     self.tol_mode, self.tol_floor, self.strength_mode,
                     spline_type='group')
        if len(self.graph_zone.page_dict[len(self.graph_zone.page_dict)]) == 9:
            self.graph_zone.page_dict[len(self.graph_zone.page_dict) + 1] = []
            self.graph_zone.num_pages += 1
            self.graph_zone.page_total.configure(text='/ %s'
                                                 % self.graph_zone.num_pages)
        self.graph_zone.page_dict[len(self.graph_zone.page_dict)].append(
            len(self.individual_dict))
        self.graph_zone.deselect_mini_graph()
        self.graph_zone.current_slot = ''
        self.current_page.set(len(self.graph_zone.page_dict))
        self.graph_zone.mini_graphs(len(self.graph_zone.page_dict))

    def open_smoothing_file(self, event=None):
        file_opt = options = {}
        options['defaultextension'] = '.csv'
        options['filetypes'] = [('all files', '.*'), ('csv files', '.csv'),
                                ('text files', '.txt')]
        options['parent'] = self.root
        options['title'] = 'Select a file...'
        spfile = filedialog.askopenfile(mode='r', **file_opt)
        if spfile is not None:
            for k in self.sp_dict.keys():
                self.sp_dict[k].set('-1')
            lines = spfile.readlines()
            spfile.close()
            for l in lines:
                tempind = int(l.split(',')[0])
                if tempind in self.individual_dict.keys():
                    newsp = str(l.split(',')[1][:-1])
                    self.sp_dict[tempind].set(newsp)
                    self.individual_dict[tempind].sp_status = 'cyan'
                    self.individual_dict[tempind].update()
            if self.graph_zone.view == 'mini':
                self.graph_zone.mini_graphs(self.current_page.get(),
                                            and_deselect=False)
                self.graph_zone.fig.canvas.draw()
            elif self.graph_zone.view == 'mega':
                self.graph_zone.update_mega_graph()
                self.update_summary()
                self.current_sp.set(self.individual_dict[
                    self.current_col.get()].smoothing_value)

    def save_smoothing_values(self, event=None):
        if platform == 'win32':
            ext = ''
        else:
            ext = '.csv'
        spfile = filedialog.asksaveasfile(mode='w',
                                          initialfile='smoothing.csv',
                                          defaultextension=ext,
                                          filetypes=[('all files', '.*'),
                                                     ('csv files', '.csv')],
                                          parent=self.root,
                                          title='Save smoothing values')
        if spfile is not None:
            for i in self.individual_dict.values():
                if i.sp_status == 'cyan':
                    spfile.write('%d,%s\n'
                                 % (i.id_number, i.smoothing_value.get()))
            spfile.close()

    def clear_smoothing_values(self, event=None):
        for i in self.individual_dict:
            if self.individual_dict[i].sp_status == 'cyan':
                temp_individual_id = self.individual_dict[i].id_number
                self.sp_dict[temp_individual_id].set('-1')
                self.individual_dict[i].update()
                self.individual_dict[i].sp_status = 'magenta'
        self.graph_zone.mini_graphs(self.current_page.get(),
                                    and_deselect=False)
        self.graph_zone.fig.canvas.draw()

    def load_settings(self, event=None):
        usrSett = shelve.open('UserSettings')
        if len(usrSett) > 0:
            self.view_names.set(usrSett['view_names'])
            self.view_pts.set(usrSett['view_pts'])
            self.view_pandtol.set(usrSett['view_pandtol'])
            self.view_spline.set(usrSett['view_spline'])
            self.view_se.set(usrSett['view_se'])
            self.tol_type.set(usrSett['tol_type'])
            self.tol_drop.set(usrSett['tol_drop'])
            self.tol_absolute.set(usrSett['tol_absolute'])
            self.tol_mode.set(usrSett['tol_mode'])
            self.tol_floor.set(usrSett['tol_floor'])
            self.loc_peak.set(usrSett['loc_peak'])
            self.peak_min.set(usrSett['peak_min'])
            self.peak_max.set(usrSett['peak_max'])
            self.strength_mode.set(usrSett['strength_mode'])
            self.sp_lim.set(usrSett['sp_lim'])
            self.sp_min.set(usrSett['sp_min'])
            self.sp_max.set(usrSett['sp_max'])
        usrSett.close()
        self.control_panel.smoothing_limits_box.sp_lim_toggle(andupdate=False)
        self.control_panel.peak_box.loc_peak_toggle(andupdate=False)
        self.control_panel.tolerance_box.change_tol_type(andupdate=False)
        self.update_magenta_graphs()

    def save_settings(self, event=None):
        usrSett = shelve.open('UserSettings')
        usrSett['view_names'] = self.view_names.get()
        usrSett['view_pts'] = self.view_pts.get()
        usrSett['view_pandtol'] = self.view_pandtol.get()
        usrSett['view_spline'] = self.view_spline.get()
        usrSett['view_se'] = self.view_se.get()
        usrSett['tol_type'] = self.tol_type.get()
        usrSett['tol_drop'] = self.tol_drop.get()
        usrSett['tol_absolute'] = self.tol_absolute.get()
        usrSett['tol_mode'] = self.tol_mode.get()
        usrSett['tol_floor'] = self.tol_floor.get()
        usrSett['loc_peak'] = self.loc_peak.get()
        usrSett['peak_min'] = self.peak_min.get()
        usrSett['peak_max'] = self.peak_max.get()
        usrSett['strength_mode'] = self.strength_mode.get()
        usrSett['sp_lim'] = self.sp_lim.get()
        usrSett['sp_min'] = self.sp_min.get()
        usrSett['sp_max'] = self.sp_max.get()
        usrSett.close()

    def reset_settings(self, event=None):
        self.settings_to_default()
        self.control_panel.smoothing_limits_box.sp_lim_toggle(andupdate=False)
        self.control_panel.peak_box.loc_peak_toggle(andupdate=False)
        self.control_panel.tolerance_box.change_tol_type(andupdate=False)
        self.update_magenta_graphs()

    def output_graphs(self, event=None):
        '''Create a pdf, svg or eps file via R of all the graphs.
        This function pays attention to the current View settengs, and so if
        data points are toggled off in the PFunc GUI, they will be absent from
        this output as well.
        '''
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        if platform == 'win32':
            ext = ''
        else:
            ext = '.pdf'
        graphfile = filedialog.asksaveasfile(mode='w',
                                             initialfile='spline_graphs.pdf',
                                             defaultextension=ext,
                                             filetypes=[('all files', '.*'),
                                                        ('pdf files', '.pdf'),
                                                        ('eps files', '.eps'),
                                                        ('svg files', '.svg')],
                                             parent=self.root,
                                             title='Select a file...')
        if graphfile is not None:
            if graphfile.name[-4:] == '.svg':
                filetype_for_r = 'svg'
            elif graphfile.name[-4:] == '.eps':
                filetype_for_r = 'eps'
            else:
                filetype_for_r = 'pdf'
            if filetype_for_r == 'pdf':
                r('''pdf(file = '%s', onefile = TRUE)
                     par(mfrow = c(3, 3), mar = c(1.5, 1.1, 2, 1.1),
                         oma = c(1, 1.5, 0, 0.5))
                     min.resp <- %s
                     max.resp <- %s
                     resp.range <- max.resp - min.resp
                     max.y <- max.resp + 0.02 * resp.range
                ''' % (graphfile.name, self.individual_dict[1].axes_ranges[2],
                       self.individual_dict[1].axes_ranges[3]))
                # isn't there a better way to handle min and max resp?
            elif filetype_for_r == 'svg':
                nrows = ceiling(len(self.individual_dict) / 3)
                svg_height = str(nrows * (7/3))
                r('''svg(file = '%s', height = %s)
                     par(mfrow = c(%s, 3), mar = c(1.5, 1.1, 2, 1.1),
                         oma = c(1, 1.5, 0, 0.5))
                         min.resp <- %s
                         max.resp <- %s
                         resp.range <- max.resp - min.resp
                         max.y <- max.resp + 0.02 * resp.range
                    ''' % (graphfile.name, svg_height, nrows,
                           self.individual_dict[2].axes_ranges[2],
                           self.individual_dict[1].axes_ranges[3]))
            elif filetype_for_r == 'eps':
                nrows = ceiling(len(self.individual_dict) / 3)
                eps_height = str(nrows * (7/3))
                r('''setEPS()
                     postscript(file = '%s', height = %s, width = 7,
                                paper = 'special')
                     par(mfrow = c(%s, 3), mar = c(1.5, 1.1, 2, 1.1),
                         oma = c(1, 1.5, 0, 0.5))
                         min.resp <- %s
                         max.resp <- %s
                         resp.range <- max.resp - min.resp
                         max.y <- max.resp + 0.02 * resp.range
                    ''' % (graphfile.name, eps_height, nrows,
                           self.individual_dict[2].axes_ranges[2],
                           self.individual_dict[1].axes_ranges[3]))
            for i in self.individual_dict:
                tempind = self.individual_dict[i]
                self.draw_one_graph_in_r(self.individual_dict[i])
            r('dev.off()')
            graphfile.close()
            self.root.config(cursor='')

    def draw_one_graph_in_r(self, individual):
        individual.update()
        isSubmerged = individual.tolerance_height > individual.peak_resp
        if self.tol_mode.get() == 'broad':
            current_tolerance_points = (
                individual.broad_tolerance_points.r_repr())
        elif self.tol_mode.get() == 'strict':
            current_tolerance_points = (
                individual.strict_tolerance_points.r_repr())
        r('''individual_data <- %s
             peak_bundle <- list(peak.response = %s,
                                 peak.preference = %s,
                                 predicting.stimuli = data.frame(stim=%s),
                                 predicted.response = %s,
                                 predicted.se = %s)
             tolerance_bundle <- list(tolerance.height = %s,
                                      cross.points = %s,
                                      submerged = %s)
             ghost_bundle <- list()
             is.flat <- CheckForFlat(individual_data, 2)
             #is.flat <- CheckForFlat(#s, 2)
             #if (sd(#s) == 0) {flat <- TRUE}
        ''' % (individual.r_data_frame.r_repr(),
               individual.peak_resp,
               individual.peak_pref,
               individual.spline_x.r_repr(),
               individual.spline_y.r_repr(),
               individual.se.r_repr(),
               individual.tolerance_height,
               current_tolerance_points,
               #individual.tolerance_points.r_repr(),
               str(isSubmerged).upper(),
               #individual.data_y.r_repr()
               ))
        if self.view_names.get() == 1:
            name = individual.name
        else:
            name = ''
        r("""plot(NA, NA, main = "%s", xlab = "", ylab = "",
                  ylim = c(%s, %s), xlim = c(%s, %s), type = "l")
          """
          % (name,
             individual.axes_ranges[2], individual.axes_ranges[3],
             individual.axes_ranges[0], individual.axes_ranges[1]))
        groupcheck = 'FALSE'
        if individual.type == 'group':
            groupcheck = 'TRUE'
        r('''GraphSpline(individual_data, peak_bundle, tolerance_bundle,
                         '%s', 2, %s,
                         %s, %s, TRUE, '%s', max.y,
                         FALSE, ghost_bundle, is.flat, %s,
                         2, forgui = TRUE, group = %s, graph.se = %s,
                         graph.spline = %s)
        ''' % (name, self.view_pts.get(),
               self.view_pandtol.get(), self.view_pandtol.get(),
               self.tol_mode.get(), individual.smoothing_value.get(),
               groupcheck, self.view_se.get(), self.view_spline.get()))

    def output_summaries(self, event=None):
        '''Output a csv file with all of the spline measures listed in the
        Summary box (peak preference, peak height, tolerance, etc.) for all
        individuals in the dataset.
        '''
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        if platform == 'win32':
            ext = ''
        else:
            ext = '.csv'
        summfile = filedialog.asksaveasfile(mode='w',
                                            initialfile='spline_summaries.csv',
                                            defaultextension=ext,
                                            filetypes=[('all files', '.*'),
                                                       ('csv files', '.csv')],
                                            parent=self.root,
                                            title='Save spline summaries...')
        if summfile is not None:
            r('''output <- data.frame(name = rep(NA, %s),
                           peak_pref=NA, peak_height=NA, tolerance=NA,
                           strength=NA,
                           #HD_strength=NA, HI_strength=NA,
                           responsiveness=NA, smoothing=NA)
            ''' % len(self.individual_dict))
            for i in self.individual_dict:
                tempind = self.individual_dict[i]
                tempind.update()
                if self.tol_mode.get() == 'broad':
                    temp_tolerance = self.individual_dict[i].broad_tolerance
                elif self.tol_mode.get() == 'strict':
                    temp_tolerance = self.individual_dict[i].strict_tolerance
                if self.strength_mode.get() == 'Height-Dependent':
                    temp_strength = tempind.hd_strength
                elif self.strength_mode.get() == 'Height-Independent':
                    temp_strength = tempind.hi_strength
                r('''output[%s, 1] <- '%s'
                     output[%s, 2:7] <- c(%s, %s, %s, %s, %s, %s)
                ''' % (i, tempind.name, i, tempind.peak_pref,
                       tempind.peak_resp,
                       temp_tolerance, temp_strength,
                       tempind.responsiveness, tempind.smoothing_value.get()))
            r("write.csv(output, '%s', row.names = FALSE)" % summfile.name)
            summfile.close()
            self.root.config(cursor='')

    def output_points(self, event=None):
        '''Output a csv file of points that make up the splines in every graph.
        The continuous curves of the splines are broken into 200 equally spaced
        points. These points can then be used to plot splines in other
        programs. x- and y-values are output for each individual, and if the
        Standard Error setting is toggled on in the View settings, then
        standard error points of the spline are output as well.
        '''
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        if platform == 'win32':
            ext = ''
        else:
            ext = '.csv'
        pointfile = filedialog.asksaveasfile(mode='w',
                                             initialfile='spline_points.csv',
                                             defaultextension=ext,
                                             filetypes=[('all files', '.*'),
                                                        ('csv files', '.csv')],
                                             parent=self.root,
                                             title='Select a file...')
        if pointfile is not None:
            r("output <- data.frame(x = rep(NA, 201))")
            for i in self.individual_dict:
                tempind = self.individual_dict[i]
                tempind.update()
                r('''output$%s_stimulus <- %s
                     output$%s_response <- %s''' % (tempind.name,
                                                    tempind.spline_x.r_repr(),
                                                    tempind.name,
                                                    tempind.spline_y.r_repr()))
                if self.view_se.get() == 1:
                    r('output$%s_se <- %s'
                        % (tempind.name, tempind.se.r_repr()))
            r('output <- output[2:ncol(output)]')
            r('write.csv(output, "%s", row.names = FALSE)' % pointfile.name)
            pointfile.close()
            self.root.config(cursor='')

    def output_tol(self, event=None):
        '''Tolerance is the width of the spline at a certain height. In the
        graphs, it is represented by a horizontal blue line. Tolerance points
        are the start and stop points of those blue lines. This function
        outputs a csv file of these tolerance points for each individual.
        Like in the output_points function, this is useful for plotting splines
        in another program.
        '''
        try:
            self.root.config(cursor='wait')
        except:
            self.root.config(cursor='watch')
        if platform == 'win32':
            ext = ''
        else:
            ext = '.csv'
        pointfile = filedialog.asksaveasfile(
            mode='w', initialfile='tolerance_points.csv', defaultextension=ext,
            filetypes=[('all files', '.*'), ('csv files', '.csv')],
            parent=self.root, title='Select a file...')
        if pointfile is not None:
            output_tol_table = ''
            for i in range(1, len(self.individual_dict) + 1):
                individual_name = self.individual_dict[i].name
                if self.tol_mode.get() == 'broad':
                    individual_tol_pts = (
                        self.individual_dict[i].broad_tolerance_points)
                elif self.tol_mode.get() == 'strict':
                    individual_tol_pts = (
                        self.individual_dict[i].strict_tolerance_points)
                tol_pts_str = ''
                for i in individual_tol_pts:
                    tol_pts_str = tol_pts_str + str(i) + ', '
                tol_pts_str = tol_pts_str[: -2]
                output_row = individual_name + ', ' + tol_pts_str
                if i != (len(self.individual_dict) + 1):
                    output_row += '\n'
                output_tol_table += output_row
            pointfile.write(output_tol_table)
            pointfile.close()
            self.root.config(cursor='')

    def quit(self, event=None):
        self.root.quit()

    def create_about_window(self, event=None):
        self.about_window = AboutWindow(self, self.about_font1)

if __name__ == '__main__':
    main_app = MainApp()
    main_app.root.mainloop()


# CHANGE LOG
#
# 151105 - Completed first basic working version after ~36 hours
# 151108 - Implemented the [open file] dialog, got it working with data
#        - Implemented menubar
# 151109 - peak.between is no longer on by default
#        - sp.binding is now on by default, between 0.05 and 5
#        - Installed set and reset buttons for sp
#        - Full functionality for sp set button
#        - Added quit button to menu bar
#        - Open and use sp files
# 151110 - Automatically "Spline it!" after loading data or sp files.
#        - Automatically clear the sp file when loading a new data file.
#        - Removed quit from File menu
#        - Full functionality of sp reset button
#        - Assign and save smoothing parameters
#        - Implemented clear all sps in File menu
#        - Added welcome screen
#        - Removed automatically loading data file
#        - Disable certain items on startup. Enable them upon opening data.
#        - Export pdf of all splines
#        - Export summary file for all splines
#        - Export points along the spline and standard error
# 151111 - Made it okay to cancel in the file dialog
#        - Automatically allow for read.csv OR read.delim
#        - Added error message for data that have less than 2 columns
#        - Program now remembers dir path for opened files
#        - Swapped position of + and - buttons
# 151117 - Keypresses only have effect when certain windows have focus
#        - Tolerance, peak within and sp limit options now functional
# 151118 - Implemented View menu
#        - Added option to not view data points
# 151124 - Lots of work over the past few days to migrate from graphing with an
#          X11 popup window from R to a 100% Python interface using matplotlib
#        - With this switch comes implemention of OOP where appropriate
#        - Graphing window is now connected to control widgets
#        - 3x3 view of graphs now possible
# 151125 - Implemented ability to select mini graphs and display stats
# 151129 - Added page nav buttons below mini graph view
# 151130 - Added full functionality to page nav buttons
#        - Individual's name and sp appear in appropriate boxes when selected
#        - Fixed bug related to re-selecting a recently deselected individual
#        - Smoothing parameters can now be adjusted in mini-graph view
#        - Changing smoothing parameters is stickier now-- don't need to press
#          "set" button to apply them.
# 151131 - Implemented 'mega' graph view-- double-click on a mini graph to
#          view it large. Smoothing adjustments work with it.
# 151201 - Made control panel an object class
#        - Started using LabelFrame in control panel
#        - Got rid of "set" sp button
# 151202 - Save, Open and Clear sp all work now.
#        - Pressing Return in the sp box sets sp.
#        - Moved summary zone to bottom of side pane
#        - Display column name in lower right graph_zone in italics
#        - Page controls now double as individual nav controls in mega view
#        - When toggling View Points, the change is immediately visualized
# 151203 - Added options to toggle splines, peak lines, and tolerance lines.
#        - Consolidated draw_graph function
# 151206 - Implemented full ability to read in vertical or horizontal data
#        - Added Tolerance controls to sidebar
# 151210 - Added SP Limits and Peak options to sidebar
#        - Disable File menu options until relevant
#        - Option to view individual IDs over graphs
#        - Removed individual name from lower right corner
# 151224 - Removed Options button.
# 160409 - Removed Quit button from menu bar (still exists uner File)
#        - Removed Ghost and SE options from View (commented out)
# 160412 - Changed root.title to PFunc
#        - Finished implementing save/load settings using shelve
#        - Fixed: re-highlight selected graph after changing settings
#        - Polished the spline points output function
#        - Linked up to new version of R script (160412m)
# 160413 - Cleaned up code
#        - Fixed bug that prevented inputting custom SPs with Return key
#        - Rearranged import statements so that they all work with each other
#        - All three output options work now... in horizontal mode
# 160414 - All three output options now work in vertical mode too
# 160428 - GUI can now handle flat data (all responses the same value)
#        - Fixed bug that prevented nav buttons from returning to first figure
#        - Fixed bug that kept old slot selected after opening new file
#        - Widened the page entry box (page_num_ent) from 2 to 3
#        - Fixed bug with using the page entry box for graph nav in mega view
#        - Added platform-specific formatting
# 160621 - Implemented tolerance floor option
#        - Red line for peak now goes to x-axis, not to 0
# 160622 - Fixed bug that kept settings from carrying over to other pages. I
#          did this by re-splining all 9 individuals every time a page turns.
#          There is definitely room for optimization here (e.g. first check
#          to see whether any settings updated before bothering to respline)
#        - Added a method to the PrefFunc class, update()
# 160623 - Fixed bugs associated with changing sp and settings mega view after
#          having navigated
#        - Added indicator lights to graphs, showing whether the sp has been
#          changed
#        - Changing smoothing parameter only updates selected graph
#        - Changes in sp limits only affects magenta graphs
#        - Clicking sp limit checkbox forces update
#        - Clicking tolerance mode radiobuttons forces update
#        - Fixed all three output functions so that they match sp and settings
# 160624 - Added message log to control panel
#        - Added smoothing to the summary window so that users can confirm
#          changes
# 160627 - Refined message log
#        - Unified horizontal and vertical output of pdf (consider doing for
#          other outputs too)
# 160628 - Created interface for higher-order splines
# 160630 - Created backend for higher-order splining
# 160701 - Addressed several bugs in higher-order splining
# 160702 - Fixed bug preventing summary box from updating after resetting sp
#        - Consolidated output summary function to be the same for horizontal
#          as for vertical. Now works with group splines.
# 160706 - Output graphs now include higher-order splines
#        - Output points now include higher-order splines (now same function
#          for vertical and horizontal)
# 160707 - Can now successfully process higher-order splines that don't
#          perfectly match up (still need to weigh the GAMs)
# 160713 - Implemented the view stanard error option in the interface
#        - Implemented the standard error option for graph and point outputs
#        - Fixed a small bug related to group-level splining
#        - Fixed a small bug related to outputting points
#        - Added a function to higher-level splining that replaces spaces with
#          underscores if they appear in the name.
# 160819 - Added "default" label to the mean option in group-level splining
#        - Changed "higher-order" splining to "group-level"
# 160824 - Can now handle missing values
#        - Can now handle data with fewer than 10 points
# 161010 - Figured out how to run PFunc in Windows. Made slight adjustment so
#          that Windows can make use of rpy2 (added a conditional statement
#          that changes "\\" in paths to "/" -- needed for passing setwd to R)
# 161011 - Added startup formatting conditionals for Windows
# 161012 - Figured out how to give PFunc a proper icon!
#        - Edited the Tolerance box on the settings panel to make room for
#          the implementation of absolute tolerance.
# 161017 - Continued to reformat the Tolerance settings box for the absolute
#          option.
#        - Fully implemented the absolute tolerance (i.e., tolerance at set
#          y-axis value). Note: it is simply a matter of setting drop to 1 and
#          floor to your desired value.
#        - Added update_summary as a method to the control panel object. Must
#          continue to migrate to using this method exclusively.
#        - Changed "Peak + Tolerance" to "Peak & Tolerance"
#        - Added get_name and get_tol_pts_str methods to PrefFunc class
#        - Added ability to output the tolerance points--that is, the x-axis
#          values corresponding to the boundaries of tolerance.
# 161019 - Added ellipsis to group spline menu option.
# 161024 - Improved formatting on Windows
#        - Improved formatting on Mac OS X
#        - Added feature to specify the location of R on the system in case
#          rpy2 cannot find it (this is in the PFuncPath.txt file).
# 161028 - Minor improvements to interface on Mac and Linux
#        - Fixed bug that appended file extension twice in Windows (although
#          if Windows users do not supply extension, file will be output with
#          no extension at all).
# 161108 - Fixed bug that prevented PFunc from running when user chose to
#          "open with" Python in Windows.
# 161116 - Minor text fixes
# 161122 - Changed the summary window: now only one measure of strength.
#        - Added setting to change how strength is measured (HI vs HD)
# 161129 - Removed Message box from side bar.
#        - Changed "SP Limits" to "Smoothing Limits"
#        - Changed order of Peak Tolerance and Strength on Settings panel
#        - Changed peak settings input to specify x-axis range, rather than
#          central proportion of range, which was less intuitive and more
#          restrictive.
# 161130 - Changed peak_win to peak_btwn1 and peak_btwn2
#        - Fully implemented the new interface for finding local peaks
#        - Increased font size of headings slightly
#        - Increased font size of summary window slightly for Windows
# 161202 - Fully implemented ability to switch between HD and HI strength
#        - Graphs resize with windows now!
#        - Decided to stick with grid over pack for geometry manager
# 161205 - Cleaned up some bugs with resizing graphs and window
#        - Added copyright statement to start page
#        - Fixed problem in Windows where data definer for vertical data wasn't
#          appearing in the middle of the window.
#        - Added About window to Help menu
#        - Added option in Help menu to open README.txt in the system's default
#          text editor
# 161206 - Removed "default" from next to "mean" in group-level tool
#        - Added message log window available from the Advanced menu
# 161207 - Finished creating and debugging the new message log
#        - Added messages to show up in the message log
#        - Fixed bug associated with graph navigation. Bug was caused by the
#          fact that PFunc keeps track of the last graph that was selected in
#          order to make double-clicking smoother. To fix the problem, I made
#          PFunc reset recent_col and recent_slot when navigating in mega mode.
# 161212 - Appropriate fields are enabled or disabled, when loading or
#          resetting settings, depending on the checkboxes that are clicked
#          (smoothing limits, local peak, tolerance type)
#        - Added strength mode and tolerance floor to saved settings
#        - Changed the default name of saved smoothing parameter file from
#          newsps to smoothing, which is more user-friendly
#        - Fixed a bug that wasn't properly saving smoothing-value files
#        - Fixed a bug that gave saved smoothing-value file double extension
#          in windows
#        - Added minor feature: when new smoothing values are loaded, the
#          the affected graphs turn cyan and update. The rest remain untouched
#        - Finalized new icon and applied it to root and all Toplevel windows
#        - The group-level spline window now pops up in front of root window
#        - Changed "Local Peak" to "Find Local Peak" in the control panel
# 161213 - Minor platform-specific adjustments to code for compatibility
#        - Found a way to convert .md files to .pdf through Atom, so now the
#          help menu will bring up a pdf instead of a text file.
# 170102 - Finished conforming to pep8 style standards (with a minor exception
#          in the import statements, which cannot be helped).
# 170111 - Changed spline graph output so that the pdf shows 3 rows, not 4, to
#          match the view in PFunc.
#        - Reformatted one of the error messages.
# 170112 - Removed get_name() method from PrefFunc class
#        - Removed get_tol_pts_str() method from PrefFunc class; implemented it
#          directly in the output_tol function.
#        - Updated copyright date.
# 170113 - Removed the RowColCoords class and the gridic placement of widgets
#          in the control panel.
#        - Re-grouped all class definitions at top of file.
# 170116 - Started major overaul of code. I'm working toward better
#          better encapsulation. Today I started by making each component of
#          the control panel its own class, and I've begun shifting toward
#          hooking up all the components by calling customized virtual events
#          when certain buttons are pushed, rather than calling instances
#          directly in the code.
# 170117 - Continued progress encapsulating the program. Finished first draft
#          of all the pieces of the control panel.
# 170118 - Continued progress on the encapsulation re-write. Turned MainApp
#          into a class, rather than having the program free-floating in the
#          file.
#        - Started a clean file to sort out useful code from junk.
# 170119 - Continued progress on the encapsulation re-write.
# 170120 - Continued progress on the encapsulation re-write.
#        - Renamed peak_btwn1 and _2 to peak_min and _max
# 170125 - Continued progress on the encapsulation re-write.
#        - Got a lot of the basics up and running, like opening files and
#          adjusting smoothing parameters.
# 170126 - Continued progress on the encapsulation re-write.
#        - Changed how smoothing parameters are assigned and called.
#        - Got reset_sp working.
#        - Connected self.current_page in root to graph_zone.
#        - Hooked up view buttons to update graphs when clicked.
#        - Bound Return key to page numbers, smoothing value
#        - Discovered and fixed bug: when slot is selected and then a new
#          page number is entered, if the page doesn't have a graph in
#          the slot, PFunc ceases to function until an acceptable page is
#          returned to. I fixed this by clearing the slot when changing
#          pages.
#        - Added upper and lower limits to enter_page_num
# 170131 - Continued progress refactoring the code.
#        - Got the Smoothing Limits, Local Peak and Tolerance boxes working.
# 170201 - Continued refactoring code.
#        - Fixed bug that was affecting how graphs are selected. Related to
#          the new and_deselect argument in select_graph method
#        - Got menu options activated when new file is opened
# 170202 - Continued refactoring the code.
#        - Finished implementing open_vertical_file
#        - Got message log working. Created a dictionary of messages and codes
# 170203 - Got the combine group spline stuff working, along with the
#          loading, saving and clearing of smoothing parameters
#        - Made good progress on the load/save/clear settings options
# 170209 - Finished all the output options in the File menu
#        - Completed About window
#        - Confirmed that I have all the pieces from the previous version
# 170210 - Slight modification to how the Help menu is created
# 170217 - Fixed bug where smoothing parameter was not updating in mega view
#        - Fixed bug that came up when opening the group-level window
#        - Discovered that the gam function can handle multiple x-values per y.
#          this could change the way we do group-level splines. I added a
#          temporary option to run group splines without first taking the mean
#          or the median to see how those splines compare.
# 170222 - Corrected all the Output Spline Figures problems
#        - Created option to output spline figures as one big .svg file
# 170224 - Created messages and errors associated with missing stimuli or
#          having too few responses.
# 170309 - Group splines are now only fit directly to constituent splines
#        - Expanded head room on graphs
#        - Fixed display of group names
#        - Added restrictions to group names--any weird symbols are replaced by
#          underscores
# 170310 - Added EPS graph output option
# 170316 - Formatted PFunc for Mac (especially fonts)
#        - Established positioning rules for the messages window
#        - Removed unnecessary R commands that assigned newsps variable
# 170426 - Changed Peak Stimulus to Peak Preference in the GUI. Note that peak
#          stim is still used in the code.
#        - Added some comments to the code.
# 170503 - Added a line up by the imports to check whether matplotlib 2+ is
#          installed (rather than 1.5), and if it is, then make it look like
#          the old version. In the fullness of time, I can update the look
#          of the UI, but this is a good solution for now.
#        - Allowed for apostrophes in filenames.
#        - Ran code by pep8. Fixed everything except for E402 - placement
#          of import statements, because some commands must come before
#          further imports.
# 170508 - Changed all references to peak stimulus to peak preference
# 170511 - Changed spstatus to sp_status
#        - Fixed bug that didn't carry over smoothing limit settings to other
#          pages
# 170512 - Changed the version convention from date of release (YYMMDD) to
#          version numbers: major.minor.patch (semantic versioning),
#          starting with 0.9.0
# 170515 - Implemented update_peak and update_tolerance methods of PrefFunc
#          class to increase efficiency.
# 170516 - Fixed a bug--now tolerance updates when it's in strict mode and when
#          a peak-update results in a new peak.
#        - Made the clear_smoothing_values function more efficient--now it only
#          clears smoothing values for cyan individuals.
#        - Added a loading screen for when opening data files.
#        - Fixed a bug related to entering smoothing parameters.
#        - Added cursor changes (busy vs ready)
# 170518 - Fixed cursor bug in Linux
#        - Added warning message when PFunc can't find README.pdf
