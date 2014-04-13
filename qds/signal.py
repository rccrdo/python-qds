# Copyright (c) 2009 Riccardo Lucchese, riccardo.lucchese at gmail.com
#
# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.
#
# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:
#
#    1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
#
#    2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
#
#    3. This notice may not be removed or altered from any source
#    distribution.

class Signal():
    def __init__(self, name="unnamed signal"):
        self._name = name
        self._timeline = []
        self._values = []
        self._values_rows = None
        self._values_cols = None
    
    def __getitem__(self, at):
        return (self._timeline[at], self._values[at])

    def __len__(self):
        return len(self._timeline)

    def shape(self):
        assert len(self._timeline)
        return self._values[0].shape

    def append(self, time, value):
        # accept only 2D numpy arrays
        assert value.ndim == 2
        if not len(self._timeline):
            self._values_rows = value.shape[0]
            self._values_cols = value.shape[1]
        else:
            # force all values to have the same shape
            assert self._values_rows == value.shape[0] and \
                   self._values_cols == value.shape[1]

        self._timeline.append(time)
        self._values.append(value)

    def timeline(self):
        return self._timeline

    def upper_triang_trajectories(self):
        if not len(self._timeline):
            return []

        assert self._values_rows == self._values_cols

        trajectories = []
        n = self._values_rows
        num_trajs = (n*(n+1))/2
        for i in range(0,num_trajs):
            trajectories.append([])

        for v in self._values:
            i = 0
            for row in range(0, n):
                for col in range(row, n):
                    trajectories[i].append(v[row][col])
                    i = i+1

        return trajectories
        

