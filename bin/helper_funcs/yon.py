def y_or_n(prompt:str)->bool:
    while(True):
        usr_input = (input(prompt)).lower()

        if   usr_input in ("y", "ye", "yes"): return True
        elif usr_input in ("n", "no"):        return False
        else: print("Invalid input. (y/n)")