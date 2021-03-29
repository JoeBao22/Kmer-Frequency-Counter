import multiprocessing
import os
import random
import time

def f(q, name, value):
    print("start", (name, value, os.getpid()))
    q.put(name, value)
    time.sleep(random.randint(1, 5))
    print("end", (name, value, os.getpid()))
    
if __name__ == "__main__":
    queue = multiprocessing.Queue()
    all_process = []
    for i in range(10):
        process = multiprocessing.Process(target=f, args=(queue, i, 0))
        process.start()
        all_process.append(process)
    for p in all_process:
        p.join()
    # print(all_process)
        