subroutine get_blocksinfo(jb,je,ib,ie)
        use blocks
        use domain

        implicit none
        type(block) :: this_block
        integer :: ib,ie,jb,je
        
        this_block = get_block(blocks_clinic(1),1)
        jb = this_block%jb
        je = this_block%je
        ib = this_block%ib
        ie = this_block%ie
end subroutine get_blocksinfo
