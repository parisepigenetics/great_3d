#! /usr/bin/env python3

class Obj(object):
  __slots__ = ('i', 'l')
  def __init__(self, i):
    self.i = i
    self.l = []
all = {}
for i in range(1000000):
  all[i] = Obj(i)
