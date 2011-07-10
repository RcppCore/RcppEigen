    EIGEN_STRONG_INLINE const Scalar *begin() const
    { return m_storage.data(); }

    EIGEN_STRONG_INLINE Scalar *begin()
    { return m_storage.data(); }

    EIGEN_STRONG_INLINE const Scalar *end() const
    { return m_storage.data() + base().size(); }

    EIGEN_STRONG_INLINE Scalar *end()
    { return m_storage.data() + base().size(); }
