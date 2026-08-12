"""
Date: 06/01/2026
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
import ctypes
 
def GUI_gen():
   global out_dict
   out_dict ={}
   dark_mode = True
   try:
      ctypes.windll.shcore.SetProcessDpiAwareness(1)
   except:
      'no DPI awareness'
   start_date = pd.Timestamp("01/01/2024")
   end_date= pd.Timestamp.today()
   root = tk.Tk()
   root.title('AppalAIR Database')
   # Set up different frames for same root
   top_frame = tk.Frame(root)
   top_frame.pack(fill="x")
   if dark_mode : top_frame.configure(bg="#222425")

   mid_frame = tk.Frame(root)
   mid_frame.pack(fill="both", expand=True)
   if dark_mode : mid_frame.configure(bg="#222425")

   end_frame = tk.Frame(root)
   end_frame.pack(fill='x')
   if dark_mode : end_frame.configure(bg="#222425")

   left_mid = tk.Frame(mid_frame,width =500)
   right_mid = tk.Frame(mid_frame,width =500)

   left_mid.grid(row=0, column=0, padx=20, sticky="n")
   right_mid.grid(row=0, column=1, padx=20, sticky="n")

   mid_frame.grid_columnconfigure(0, weight=1)
   mid_frame.grid_columnconfigure(1, weight=1)
   # Create a StringVar to associate with the l1
   text_var = tk.StringVar()
   root.tk.call('tk', 'scaling', 2.0)
   root.geometry("850x1000")
   #Set up dark mode
   if dark_mode : fnt_clr = "#CEDCDE"
   else : fnt_clr = 'black'
   if dark_mode : bck_clr = "#222425"
   else : bck_clr = end_frame.cget('bg')
   if dark_mode : line_clr = "#CEDCDE"
   else : line_clr = "#333533"
   if dark_mode : srch_clr = "#CEDCDE"
   else : srch_clr = 'white'


   '''++++SET UP TITLE++++'''
   text_var.set('AppalAIR Database Access')
   l1 = tk.Label(top_frame, 
                  textvariable=text_var, 
                  anchor=tk.CENTER, justify=tk.CENTER,      
                  bg="#333533", fg="White",   
                  height=2, width=30,              
                  bd=3,                  
                  font=("Times New Roman", 16, "bold"), 
                  cursor="arrow",
                  relief=tk.GROOVE,              
                  wraplength=750         
                  )
   l1.pack(padx=15,pady=10)

   '''++++SET UP DATE RANGE SELECT++++'''
   l2 = tk.Label(top_frame, 
                  text="1: Select Date Range", 
                  anchor=tk.CENTER, justify=tk.LEFT,         
                  height=1, width=15,                            
                  font=("Times New Roman", 12, "bold"), 
                  cursor="arrow",   
                  fg=fnt_clr, bg=bck_clr,
                  wraplength=250         
                  )
   l2.pack(padx=15, pady=2)

   '''Updating date output'''
   hLeft = tk.DoubleVar(value = 0.0)  #left handle variable initialised to value 0.2
   hRight = tk.DoubleVar(value = 1.0)  #right handle variable initialised to value 0.85
   def dateConv(val, start_date,end_date):
      total_time = pd.Timedelta(end_date- start_date).total_seconds()/3600/24 #days been start and end
      date = total_time*val
      date_time = pd.Timestamp(start_date) + pd.Timedelta(date,"days")
      return pd.to_datetime(date_time).date()

   date_var = tk.StringVar(top_frame,
                           f"[ {dateConv(hLeft.get(), start_date,end_date)}                    ⟺                    {dateConv(hRight.get(), start_date,end_date)} ]")
   dates = tk.Label(top_frame, 
                  textvariable=date_var, 
                  anchor=tk.CENTER,  justify=tk.LEFT,         
                  height=0, width=50,                            
                  font=("Times New Roman", 12, "bold"), 
                  cursor="arrow",   
                  fg=fnt_clr, bg=bck_clr, 
                  padx=15, pady=0, 
                  wraplength=750         
                  )
   # dates.grid(row=2, column=0, columnspan=2)
   dates.pack(in_=top_frame)

   hSlider = RangeSliderH(top_frame, 
                        [hLeft, hRight],
                        Width= 600,
                        Height= 75, 
                        bar_color_inner = "#333533",
                        line_color = line_clr,
                        line_s_color = '#ffcc00',
                        bgColor= bck_clr,
                        digit_precision = '.3f',
                        show_value = False,
                        padX = 12)  #horizontal slider, [padX] value might be needed to be different depending on system, font and handle size. Usually [padX] = 12 serves,
                                                               #otherwise a recommended value will be shown through an error message
   hSlider.pack(in_=top_frame)   # or grid or place method could be used
   def doSomething(*args):
      date_var.set(f"{dateConv(hLeft.get(), start_date,end_date)}                ⟺                {dateConv(hRight.get(), start_date,end_date)}")

   hLeft.trace_add('write', doSomething)
   hRight.trace_add('write', doSomething)

   '''++++SET UP VARIABLE SELECT++++'''
   l3 = tk.Label(left_mid, 
                  text="2: Select Variables", 
                  anchor=tk.CENTER, justify=tk.CENTER,         
                  height=2, width=50,  
                  wraplength=500,
                  padx=15, pady=1,                             
                  font=("Times New Roman", 12, "bold"), 
                  cursor="arrow",   
                  fg=fnt_clr, bg=bck_clr,    
                  )
   #listbox options, pull this from database 
   values = ["CCN SS", "SMPS ultrafine count", "AE33 Black Carbon", "Standard Temperature","CCN flow", "SMPS fine count", "AE33 Brown Carbon", "Ambient Temperature"] 

   # create a Listbox widget
   lb_vars = tk.Listbox(left_mid, selectmode="multiple", 
                        exportselection=0, height=20)
   for value in values:
      lb_vars.insert(tk.END, value)

   #Create a scollbar
   sb_vars = tk.Scrollbar(lb_vars)
   sb_vars.pack(side=tk.RIGHT, fill=tk.Y)

   lb_vars.config(yscrollcommand=sb_vars.set) 
   sb_vars.config(command=lb_vars.yview)

   #create search functionality
   def_var = 'Search for variables...'
   search_var = tk.StringVar(value=def_var)
   vars_entry = tk.Entry(left_mid, textvariable=search_var, 
                        font=("Times New Roman", 12),fg ='grey')

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

   var_state = {"placeholder_active": True}
   def on_focus_in_var(event):
      if var_state["placeholder_active"]:
         search_var.set("")
         vars_entry.config(fg="black")
         var_state["placeholder_active"] = False
   def on_focus_out_var(event):
      if search_var.get().strip() == "":
         search_var.set(def_var)
         vars_entry.config(fg="grey")
         var_state["placeholder_active"] = True
   vars_entry.bind("<FocusIn>", on_focus_in_var)
   vars_entry.bind("<FocusOut>", on_focus_out_var)

   search_var.trace("w", update_suggestions_vars)
   lb_vars.bind("<<ListboxSelect>>", select_vars)

   l3.pack()
   vars_entry.pack()
   lb_vars.pack(side="top",padx = 10, pady = 5,expand = tk.YES, fill = "both")

   '''++++ SET UP ADDITIONAL PROCESSING SELECT ++++'''
   l4 = tk.Label(right_mid, 
                  text="3: Additional Processing", 
                  anchor=tk.CENTER,          
                  height=2,              
                  width=50,  
                  wraplength=500,
                  padx=15,               
                  pady=1,                           
                  font=("Times New Roman", 12, "bold"), 
                  cursor="arrow",   
                  fg=fnt_clr, bg=bck_clr,                
                  justify=tk.CENTER, 
                  )
   #listbox options  
   process = ["Daily Mean", 'Weekly Mean', 'Monthly Mean','STP Conversion', 'ATP Conversion'] 

   # create a Listbox widget
   lb_proc = tk.Listbox(right_mid, selectmode="multiple", 
                        exportselection=0, height=20)
   for proc in process:
      lb_proc.insert(tk.END, proc)

   #Create a scollbar
   sb_proc = tk.Scrollbar(lb_proc)
   sb_proc.pack(side=tk.RIGHT, fill=tk.Y)

   lb_proc.config(yscrollcommand=sb_proc.set) 
   sb_proc.config(command=lb_proc.yview)

   #create search functionality
   default = 'Search for processing...'
   search_proc = tk.StringVar(value = default)
   proc_entry = tk.Entry(right_mid, textvariable=search_proc,
                        font=("Times New Roman", 12), fg='grey')

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

   proc_state = {"placeholder_active": True}
   def on_focus_in_proc(event):
      if proc_state["placeholder_active"]:
         search_proc.set("")
         proc_entry.config(fg="black")
         proc_state["placeholder_active"] = False
   def on_focus_out_proc(event):
      if search_proc.get().strip() == "":
         search_proc.set(default)
         proc_entry.config(fg="grey")
         proc_state["placeholder_active"] = True
   proc_entry.bind("<FocusIn>", on_focus_in_proc)
   proc_entry.bind("<FocusOut>", on_focus_out_proc)

   search_proc.trace("w", update_suggestions_proc)
   lb_proc.bind("<<ListboxSelect>>", select_proc)

   l4.pack()
   proc_entry.pack()
   lb_proc.pack(side="top",padx = 10, pady = 5,
            expand = tk.YES, fill = "both")

   # Function for printing the
   # selected listbox value(s)
   vars_out = []
   proc_out = []
   def selected_item():
      global out_dict
      out_dict = {}
      for i in lb_vars.curselection():
         vars_out.append(lb_vars.get(i))
      for i in lb_proc.curselection():
         proc_out.append(lb_proc.get(i))
      out_dict = {'Date' : [dateConv(hLeft.get(), start_date,end_date).strftime("%d/%m/%Y"), 
                            dateConv(hRight.get(), start_date,end_date).strftime("%d/%m/%Y")],
                  'Variables' : vars_out,
                  'Processing' : proc_out}
      root.quit()

   # Create a button widget and
   # map the command parameter to
   # selected_item function
   btn = tk.Button(end_frame, text='Output Selected Data', command=selected_item)

   # Placing the button and listbox
   btn.pack(pady=5)

   lblank= tk.Label(end_frame, 
                  text="", 
                  anchor=tk.CENTER,
                  bg = bck_clr,          
                  height=2,              
                  width=15,                            
                  padx=15,               
                  pady=2,                            
                  justify=tk.LEFT,              
                  wraplength=250         
                  )
   lblank.pack(pady=15)
   root.mainloop()
   print(out_dict)
   input("end")
   return out_dict

if __name__ =='__main__':
   GUI_gen()