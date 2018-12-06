graph [
  name "island_network"
  node [
    id 0
    label "island_0"
    x_coord "30"
    y_coord "10"
  ]
  node [
    id 1
    label "island_1"
    x_coord "20"
    y_coord "30"
  ]
  node [
    id 2
    label "island_2"
    x_coord "10"
    y_coord "10"
  ]
  node [
    id 3
    label   "island_3"
    x_coord "20"
    y_coord "0"
  ]
  node [
    id 4
    label   "island_4"
    x_coord "10"
    y_coord "20"
  ]
  edge [
    id 0
    source 0
    target 1
    value 1.0
  ]
  edge [
    id 1
    source 1
    target 2
    value 0.5
  ]
  edge [
    id 2
    source 2
    target 3
    value 0.5
  ]
  edge [
    id 3
    source 4
    target 0
    value 1.0
  ]
  edge [
    id 4
    source 4
    target 3
    value 0.5
  ]
  edge
  [ id 5
    source 3
    target 1
    value 1.0
  ]
  edge
  [ id 6
    source 2
    target 4
    value 0.5
  ]
  edge [
    id 7
    source 1
    target 4
    value 1.0
  ]
]

