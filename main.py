import tkinter as tk


import MS_RIDD_GUI

# pd.options.display.width = 200

def main():
    root = tk.Tk()
    app = MS_RIDD_GUI.MSRIDD_GUI(master=root)
    app.mainloop()


if __name__ == "__main__":
    main()