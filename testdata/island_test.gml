graph [
  name "island_network"
  node [
    id 0
    label "0"
    SUID "0"
    sharedname "island_0"
    name "island_4"
    selected "false"
    x_coord "30"
    y_coord "10"
  ]
  node [
    id 1
    label "1"
    SUID "2"
    sharedname "island_1"
    name "island_1"
    selected "false"
    x_coord "20"
    y_coord "30"
  ]
  node [
    id 2
    label "2"
    sharedname "island_2"
    name "island_2"
    selected "false"
    x_coord "10"
    y_coord "10"
  ]
  node [
    id 3
    label "3"
    sharedname "island_3"
    name "island_3"
    selected "false"
    x_coord "20"
    y_coord "0"
  ]
  node [
    id 4
    label "4"
    SUID "4"
    sharedname "island_4"
    name "island_4"
    selected "false"
    x_coord "10"
    y_coord "20"
  ]
  edge [
    id 0
    source 0
    target 1
    value 1.0
    SUID "104"
    sharedname "Node 0 (interaction type) Node 1"
    sharedinteraction "interaction type"
    name "Node 0 (interaction type) Node 1"
    selected "false"
    interaction "interaction type"
  ]
  edge [
    id 1
    source 1
    target 2
    value 0.5
    SUID "103"
    sharedname "Node 1 (interaction type) Node 2"
    sharedinteraction "interaction type"
    name "Node 1 (interaction type) Node 2"
    selected "false"
    interaction "interaction type"
  ]
  edge [
    id 2
    source 2
    target 3
    value 0.5
    SUID "102"
    sharedname "Node 2 (interaction type) Node 3"
    sharedinteraction "interaction type"
    name "Node 2 (interaction type) Node 3"
    selected "false"
    interaction "interaction type"
  ]
  edge [
    id 3
    source 4
    target 0
    value 1.0
    SUID "101"
    sharedname "Node 4 (interaction type) Node 0"
    sharedinteraction "interaction type"
    name "Node 4 (interaction type) Node 0"
    selected "false"
    interaction "interaction type"
  ]
  edge [
    id 4
    source 4
    target 3
    value 0.5
    SUID "100"
    sharedname "Node 4 (interaction type) Node 3"
    sharedinteraction "interaction type"
    name "Node 4 (interaction type) Node 3"
    selected "false"
    interaction "interaction type"
  ]
  edge
  [ id 5
    source 3
    target 1
    value 1.0
    SUID "99"
    sharedname "Node 3 (interaction type) Node 1"
    sharedinteraction "interaction type"
    name "Node 3 (interaction type) Node 1"
    selected "false"
    interaction "interaction type"
  ]
  edge
  [ id 6
    source 2
    target 4
    value 0.5
    SUID "98"
    sharedname "Node 2 (interaction type) Node 4"
    sharedinteraction "interaction type"
    name "Node 2 (interaction type) Node 4"
    selected "false"
    interaction "interaction type"
  ]
  edge [
    id 7
    source 1
    target 4
    value 1.0
    SUID "97"
    sharedname "Node 1 (interaction type) Node 4"
    sharedinteraction "interaction type"
    name "Node 1 (interaction type) Node 4"
    selected "true"
    interaction "interaction type"
  ]
]
