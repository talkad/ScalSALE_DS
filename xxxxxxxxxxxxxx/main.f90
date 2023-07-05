program main
    use indexer_module

    type(indexer_t), pointer :: indexer1, indexer2

    indexer1 => get_instance(1,1,1,1)

    print*, allocated(indexer1%mapper)

    indexer1%mapper(1,0,0,0) = 50

    indexer2 => get_instance()
    print*, allocated(indexer2%mapper)

    print*, indexer1%mapper(1,0,0,0)
    print*, indexer2%mapper(1,0,0,0) 





end program main