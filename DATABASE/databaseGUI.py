"""
Date: 12/31/2025
Author: Ben Sykes
Purpose: Generate GUI for database access
"""

"""IMPORTS"""
import numpy as np
import pandas as pd 
import os
import datetime as dt 
import tkinter as tk
# from tkinter import ttk
from rangeslider import RangeSliderH 


start_date = pd.Timestamp("01/01/2024")
end_date= pd.Timestamp("01/01/2025")
root = tk.Tk()
root.title('AppalAIR Database')
# Create a StringVar to associate with the l1
text_var = tk.StringVar()
root.geometry("700x800")
text_var.set('AppalAIR Database Access')


l1 = tk.Label(root, 
                 textvariable=text_var, 
                 anchor=tk.CENTER,       
                 bg="#333533",      
                 height=2,              
                 width=30,              
                 bd=3,                  
                 font=("Times New Roman", 16, "bold"), 
                 cursor="hand2",   
                 fg="White",             
                 padx=15,               
                 pady=5,                
                 justify=tk.CENTER,    
                 relief=tk.RAISED,              
                 wraplength=250         
                )
l1.pack(pady=20)

l2 = tk.Label(root, 
                 text="1: Select Date Range", 
                 anchor=tk.CENTER,          
                 height=3,              
                 width=15,                            
                 font=("Times New Roman", 12, "bold"), 
                 cursor="hand2",   
                 fg="Black",  
                 padx=15,               
                 pady=5,                            
                 justify=tk.LEFT,              
                 wraplength=250         
                )
l2.pack(pady=5)
 
hLeft = tk.DoubleVar(value = 0.0)  #left handle variable initialised to value 0.2
hRight = tk.DoubleVar(value = 1.0)  #right handle variable initialised to value 0.85
hSlider = RangeSliderH(root, 
                       [hLeft, hRight],
                       Width= 600,
                       Height= 80, 
                       bar_color_inner = "#333533",
                       line_color = '#222222',
                       line_s_color = '#ffcc00',
                       digit_precision = '.3f',
                       show_value = False,
                       padX = 12)   #horizontal slider, [padX] value might be needed to be different depending on system, font and handle size. Usually [padX] = 12 serves,
                                                              #otherwise a recommended value will be shown through an error message
hSlider.pack()   # or grid or place method could be used
def doSomething(*args):
   # print(dateConv(hLeft.get(), start_date,end_date)) #  Print when variable changes.
   # print(dateConv(hRight.get(), start_date,end_date)) #  Print when variable changes.
   date_var.set(f"{dateConv(hLeft.get(), start_date,end_date)}                      {dateConv(hRight.get(), start_date,end_date)}")
def dateConv(val, start_date,end_date):
   total_time = pd.Timedelta(end_date- start_date).total_seconds()/3600/24 #days been start and end
   date = total_time*val
   date_time = pd.Timestamp(start_date) + pd.Timedelta(date,"days")
   return pd.to_datetime(date_time).date()

   

date_var = tk.StringVar(root,
                         f"{dateConv(hLeft.get(), start_date,end_date)}                      {dateConv(hRight.get(), start_date,end_date)}")
dates = tk.Label(root, 
                 textvariable=date_var, 
                 anchor=tk.CENTER,          
                 height=3,              
                 width=50,                            
                 font=("Times New Roman", 12, "bold"), 
                 cursor="hand2",   
                 fg="Black",  
                 padx=15,               
                 pady=5,                            
                 justify=tk.LEFT,              
                 wraplength=250         
                )
dates.pack(pady=0)

hLeft.trace_add('write', doSomething)
hRight.trace_add('write', doSomething)

  #return type list of format [ left handle value, right handle value ]
l3 = tk.Label(root, 
                 text="2: Select Variables", 
                 anchor=tk.CENTER,          
                 height=1,              
                 width=15,                            
                 font=("Times New Roman", 12, "bold"), 
                 cursor="hand2",   
                 fg="Black",  
                 padx=15,               
                 pady=2,                            
                 justify=tk.LEFT,              
                 wraplength=250         
                )
#listbox options  
values = ["CCN SS", "SMPS ultrafine count", "AE33 Black Carbon", "Standard Temperature","CCN flow", "SMPS fine count", "AE33 Brown Carbon", "Ambient Temperature"] 

# create a Listbox widget
lb_vars = tk.Listbox(root, selectmode="multiple", 
                     exportselection=0, height=2)
for value in values:
   lb_vars.insert(tk.END, value)

#Create a scollbar
sb_vars = tk.Scrollbar(lb_vars)
sb_vars.pack(side=tk.RIGHT, fill=tk.Y)

lb_vars.config(yscrollcommand=sb_vars.set) 
sb_vars.config(command=lb_vars.yview)

#create search functionality
search_var = tk.StringVar()
vars_entry = tk.Entry(root, textvariable=search_var, font=("Times New Roman", 12, "bold"))

store_slct_vars = []

def update_suggestions_vars(*args):
    
   search_term = search_var.get()
   suggestions = values

   matching_suggestions = [suggestion for suggestion in suggestions if suggestion.lower().startswith(search_term.lower())]

   lb_vars.delete(0, tk.END)
   for suggestion in matching_suggestions:
      lb_vars.insert(tk.END, suggestion)
   for suggestion in suggestions:
      if suggestion not in matching_suggestions:
         lb_vars.insert(tk.END, suggestion)
   for slct in store_slct_vars:
      i = lb_vars.get(0, "end").index(slct)
      lb_vars.select_set(i)
   

def select_vars(event):
   global store_slct_vars
   store = []
   selected = lb_vars.curselection()
   for var in selected: 
      keep = lb_vars.get(var)
      if keep not in store:
         store.append(keep)
         store_slct_vars = store


search_var.trace("w", update_suggestions_vars)
lb_vars.bind("<<ListboxSelect>>", select_vars)

l3.pack(pady=0)
vars_entry.pack()
lb_vars.pack(side="top",padx = 10, pady = 5,expand = tk.YES, fill = "both")


l4 = tk.Label(root, 
                 text="3: Additional Processing", 
                 anchor=tk.CENTER,          
                 height=1,              
                 width=15,                            
                 font=("Times New Roman", 12, "bold"), 
                 cursor="hand2",   
                 fg="Black",  
                 padx=15,               
                 pady=2,                            
                 justify=tk.LEFT,              
                 wraplength=250         
                )
#listbox options  
process = ["Daily Mean", 'Weekly Mean', 'Monthly Mean','STP Conversion', 'ATP Conversion'] 

# create a Listbox widget
lb_proc = tk.Listbox(root, selectmode="multiple", 
                     exportselection=0, height=2)
for proc in process:
   lb_proc.insert(tk.END, proc)

#Create a scollbar
sb_proc = tk.Scrollbar(lb_proc)
sb_proc.pack(side=tk.RIGHT, fill=tk.Y)

lb_proc.config(yscrollcommand=sb_proc.set) 
sb_proc.config(command=lb_proc.yview)

#create search functionality
search_proc = tk.StringVar()
proc_entry = tk.Entry(root, textvariable=search_proc, font=("Times New Roman", 12, "bold"))

store_slct_proc = []

def update_suggestions_proc(*args):
    
   search_term = search_proc.get()
   suggestions = process

   matching_suggestions = [suggestion for suggestion in suggestions if suggestion.lower().startswith(search_term.lower())]

   lb_proc.delete(0, tk.END)
   for suggestion in matching_suggestions:
      lb_proc.insert(tk.END, suggestion)
   for suggestion in suggestions:
      if suggestion not in matching_suggestions:
         lb_proc.insert(tk.END, suggestion)
   for slct in store_slct_proc:
      i = lb_proc.get(0, "end").index(slct)
      lb_proc.select_set(i)
   

def select_proc(event):
   global store_slct_proc
   store = []
   selected = lb_proc.curselection()
   for proc in selected: 
      keep = lb_proc.get(proc)
      if keep not in store:
         store.append(keep)
         store_slct_proc = store

search_proc.trace("w", update_suggestions_proc)
lb_proc.bind("<<ListboxSelect>>", select_proc)

l4.pack(pady=0)
proc_entry.pack()
lb_proc.pack(side="top",padx = 10, pady = 5,
          expand = tk.YES, fill = "both")

# Function for printing the
# selected listbox value(s)
vars_out = []
proc_out = []
def selected_item():
   for i in lb_vars.curselection():
      vars_out.append(lb_vars.get(i))
   for i in lb_proc.curselection():
      proc_out.append(lb_proc.get(i))
   print(vars_out)
   print(proc_out)
   root.quit()

# Create a button widget and
# map the command parameter to
# selected_item function
btn = tk.Button(root, text='Output Selected Data', command=selected_item)

# Placing the button and listbox
btn.pack(pady=5)

lblank= tk.Label(root, 
                 text="", 
                 anchor=tk.CENTER,          
                 height=2,              
                 width=15,                            
                 padx=15,               
                 pady=2,                            
                 justify=tk.LEFT,              
                 wraplength=250         
                )
lblank.pack(pady=15)
root.mainloop()
input("end")