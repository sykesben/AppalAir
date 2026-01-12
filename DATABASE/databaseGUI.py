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
from tkinter import ttk
from rangeslider import RangeSliderH 


start_date = pd.Timestamp()
root = tk.Tk()
root.title('AppalAIR Database')
# Create a StringVar to associate with the l1
text_var = tk.StringVar()
root.geometry("800x500")
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
   print(hLeft.get()) #  Print when variable changes.
   print(hRight.get()) #  Print when variable changes.
   dates = 
def dateConv(float, start_date):


dates = tk.Label(root, 
                 text=f"{dateConv(hLeft.get(), start_date)}-{dateConv(hRight.get(), start_date)}", 
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

l3.pack(pady=0)
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

l4.pack(pady=0)
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