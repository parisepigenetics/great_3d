#! /usr/bin/env python3

class Obj(object):
  def __init__(self, i):
    self.i = i
    self.l = []
all = {}
for i in range(1000000):
  all[i] = Obj(i)
